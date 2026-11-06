"""Natron script: linear Grade, protected labels, editable project and graph audit.

Set CELLBIO_COMPOSITE_JOB to a project JSON. Prepare with NatronRenderer -t,
then reopen scene.ntp with this script as -l (job action=audit) before rendering.
Reopening initializes ReadOIIO's layer menu in the qualified Natron 2.5 build.
Final display conversion uses the saved Blender color settings separately.
"""
import hashlib
import json
import math
import os
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def choice(node, name, value):
    parameter = node.getParam(name)
    parameter.setValue(parameter.getOptions().index(value))


def io_settings(node, config, reader):
    node.getParam('ocioConfigFile').setValue(str(config))
    for key in ('ocioInputSpace', 'ocioOutputSpace'):
        node.getParam(key).setValue('scene_linear')
    node.getParam('ocioInputSpaceSet' if reader else 'ocioOutputSpaceSet').setValue(True)
    if reader:
        choice(node, 'filePremult', 'premult'); choice(node, 'outputPremult', 'premult')
    else:
        choice(node, 'inputPremult', 'premult'); choice(node, 'formatType', 'input')
        choice(node, 'outputComponents', 'RGBA'); choice(node, 'bitDepth', '32f')


def graph_state(app):
    result = {'project': {'format': app.getProjectParam('outputFormat').getValue(),
                         'fps': app.getProjectParam('frameRate').getValue(),
                         'frames': [app.getProjectParam('frameRange').getValue(i) for i in range(2)]}}
    fields = {'BeautySource': ['filename','decodingPluginID','firstFrame','lastFrame','frameMode','timeOffset','filePremult','outputPremult','ocioConfigFile','ocioInputSpace','ocioOutputSpace'],
              'ProtectedLabels': ['filename','decodingPluginID','firstFrame','lastFrame','frameMode','timeOffset','filePremult','outputPremult','ocioConfigFile','ocioInputSpace','ocioOutputSpace'],
              'ArtisticGrade': ['premult','clampBlack','clampWhite','NatronOfxParamProcessA'],
              'LabelsOverBeauty': ['operation'],
              'LinearOutput': ['filename','encodingPluginID','inputPremult','formatType','outputComponents','bitDepth','ocioConfigFile','ocioInputSpace','ocioOutputSpace']}
    for name, parameters in fields.items():
        node = app.getNode(name)
        values = {key: node.getParam(key).getValue() for key in parameters}
        if name == 'ArtisticGrade':
            values.update({key: [node.getParam(key).getValue(i) for i in range(4)]
                           for key in ('blackPoint','whitePoint','black','white','multiply','offset','gamma')})
        result[name] = {'plugin': node.getPluginID(), 'parameters': values,
                        'inputs': {node.getInputLabel(i): node.getInput(i).getScriptName()
                                   for i in range(node.getMaxInputCount()) if node.getInput(i) is not None}}
    return result


def createInstance(app, group):
    job_path = Path(os.environ['CELLBIO_COMPOSITE_JOB'])
    job = json.loads(job_path.read_text()); output = Path(job['output'])
    handoff_path = Path(job['handoff']); handoff = json.loads(handoff_path.read_text()); inputs = handoff_path.parent
    if handoff['kind'] != 'linear-composite-handoff' or handoff['alpha'] != 'premultiplied' or handoff['colors']['linear_space'] != 'Linear Rec.709':
        raise ValueError('Require a qualified linear Rec.709 premultiplied handoff')
    for name, expected in handoff['files_sha256'].items():
        if digest(inputs/name) != expected: raise ValueError('Input changed: '+name)
    config = Path(job['ocio_config'])
    if digest(config) != job['ocio_sha256']: raise ValueError('Natron OCIO config changed')
    grade = job['grade']
    for key in ('multiply', 'gamma'):
        if len(grade[key]) != 3 or any(not math.isfinite(v) or v <= 0 for v in grade[key]):
            raise ValueError('Provide three finite positive '+key+' values')
    if job.get('action') == 'audit':
        expected = json.loads((output/'natron-graph.json').read_text())
        if graph_state(app) != expected['graph']: raise ValueError('Saved Natron graph differs from the prepared recipe')
        (output/'reopened-audit.json').write_text(json.dumps({'passed': True, 'project_sha256': digest(output/'scene.ntp'),
            'handoff_sha256': digest(handoff_path), 'graph': expected['graph']}, indent=2)+'\n')
        return
    if output.exists(): raise ValueError('Choose a fresh Natron output directory')
    output.mkdir(parents=True)
    width, height = handoff['resolution']
    format_name = job.get('project_format', 'CellBiology')
    app.addFormat('%s %d x%d 1' % (format_name, width, height))
    project_format = app.getProjectParam('outputFormat')
    if format_name not in project_format.getOptions():
        raise ValueError('Size already has a preset; set job.project_format to its ID: '+str(project_format.getOptions()))
    project_format.setValue(project_format.getOptions().index(format_name))
    readers = []
    for name, role in [('BeautySource','Beauty'), ('ProtectedLabels','Labels')]:
        reader = app.createReader(str(inputs/handoff['patterns'][role]))
        choice(reader, 'decodingPluginChoice', 'fr.inria.openfx.ReadOIIO')
        if reader.getParam('decodingPluginID').getValue() != 'fr.inria.openfx.ReadOIIO': raise ValueError('Expected the qualified EXR reader')
        reader.getParam('firstFrame').setValue(min(handoff['frames']))
        reader.getParam('lastFrame').setValue(max(handoff['frames']))
        choice(reader, 'frameMode', 'timeOffset'); reader.getParam('timeOffset').setValue(0)
        reader.setScriptName(name); io_settings(reader, config, True); readers.append(reader)
    node = app.createNode('net.sf.openfx.GradePlugin', 2); node.setScriptName('ArtisticGrade')
    for key in ('multiply', 'gamma'):
        for index, value in enumerate(grade[key]+[1.]): node.getParam(key).setValue(value, index)
    for key, value in [('premult',True),('clampBlack',False),('clampWhite',False),('NatronOfxParamProcessA',False)]:
        node.getParam(key).setValue(value)
    if not node.connectInput(0, readers[0]): raise ValueError('Cannot connect Beauty to Grade')
    merge = app.createNode('net.sf.openfx.MergePlugin', 2); merge.setScriptName('LabelsOverBeauty')
    choice(merge, 'operation', 'over')
    if not merge.connectInput(0, node) or not merge.connectInput(1, readers[1]): raise ValueError('Cannot merge protected labels')
    writer = app.createWriter(str(output/'linear_####.exr')); writer.setScriptName('LinearOutput')
    choice(writer, 'encodingPluginChoice', 'fr.inria.openfx.WriteOIIO')
    if writer.getParam('encodingPluginID').getValue() != 'fr.inria.openfx.WriteOIIO': raise ValueError('Expected the qualified EXR writer')
    io_settings(writer, config, False)
    if not writer.connectInput(0, merge): raise ValueError('Cannot connect the linear output')
    app.getProjectParam('frameRate').setValue(handoff['fps']/handoff['fps_base'])
    app.getProjectParam('frameRange').set(min(handoff['frames']), max(handoff['frames']))
    for item, point in zip([readers[0], node, readers[1], merge, writer], [(0,0),(0,120),(250,120),(0,240),(0,360)]):
        item.setPosition(*point)
    if not app.saveProjectAs(str(output/'scene.ntp')): raise ValueError('Natron project save failed')
    record = {'handoff': str(handoff_path), 'handoff_sha256': digest(handoff_path),
              'job_sha256': digest(job_path), 'ocio_sha256': digest(config), 'grade': grade,
              'frames': handoff['frames'], 'graph': graph_state(app),
              'output': 'Scene-linear EXR; apply the recorded Blender display transform once.',
              'rendered': False}
    (output/'natron-graph.json').write_text(json.dumps(record, indent=2)+'\n')
    app.closeProject()
