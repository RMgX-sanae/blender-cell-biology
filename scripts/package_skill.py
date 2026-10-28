"""Package instructions and generators; leave rendered media in project folders."""
import argparse
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=args.source.resolve();output=args.output.resolve()
    if output.is_relative_to(root):parser.error('The archive must be outside the skill directory.')
    allowed={'.md','.py','.json','.yaml','.yml','.pdb','.cif','.txt'}
    excluded={'__pycache__','frames','renders','outputs','validation','previews','.git'}
    members=[]
    for path in sorted(root.rglob('*')):
        relative=path.relative_to(root)
        if not path.is_file() or path.is_symlink():continue
        if excluded.intersection(relative.parts) or path.suffix.lower() not in allowed:continue
        if len(relative.parts)>1 and relative.parts[0] not in {'scripts','references','assets','agents'}:continue
        members.append(path)
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=output.with_suffix(output.suffix+'.tmp')
    with ZipFile(temporary,'w',ZIP_DEFLATED) as archive:
        for path in members:archive.write(path,Path(root.name)/path.relative_to(root))
    with ZipFile(temporary) as archive:
        assert archive.testzip() is None
        assert all(Path(name).suffix.lower() in allowed for name in archive.namelist())
        assert f'{root.name}/SKILL.md' in archive.namelist()
    temporary.replace(output)
    print(json.dumps({'archive':str(output),'files':len(members),'bytes':output.stat().st_size,
                      'rendered_media_included':False,'prepared_meshes_or_fonts_included':False},indent=2))


if __name__=='__main__':main()
