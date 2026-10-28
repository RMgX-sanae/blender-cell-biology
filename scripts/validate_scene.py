"""Read-only structural checks for a saved starter .blend; not a biological certification."""
import argparse
import json
import math
import sys
from pathlib import Path
import bpy
from mathutils import Vector


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report",type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index("--")+1:])
    scene=bpy.context.scene
    failures=[]
    checks=[]

    def check(condition,description):
        checks.append(description)
        if not condition:
            failures.append(description)

    check(scene.camera is not None,"camera exists")
    check(scene.frame_end>scene.frame_start,"positive animation range")
    sample_frames=sorted({1,scene.frame_end//2,scene.frame_end,round(scene.frame_end*.94)-1,round(scene.frame_end*.94)})
    for frame in sample_frames:
        scene.frame_set(frame)
        check(all(math.isfinite(v) for obj in scene.objects for row in obj.matrix_world for v in row),f"finite object transforms at frame {frame}")
        for obj in scene.objects:
            if obj.type=="MESH":
                check(all(math.isfinite(v) for vertex in obj.data.vertices for v in vertex.co),f"finite vertices: {obj.name} @ {frame}")

    receptors=[o for o in scene.objects if o.name.startswith("Receptor_")]
    if receptors:
        points=[o.location for o in receptors]
        check(min(((a-b).length for i,a in enumerate(points) for b in points[i+1:]),default=float("inf"))>=.62-1e-6,"receptor anchor exclusion distance")
        check(all(o.get("outward_normal")==[0.0,0.0,1.0] or tuple(o.get("outward_normal",[]))==(0,0,1) for o in receptors),"receptor extracellular direction +Z")

    cargo=[o for o in scene.objects if o.name.startswith("Luminal_cargo_")]
    for frame in (1,scene.frame_end):
        scene.frame_set(frame)
        for obj in cargo:
            check(obj.parent and obj.parent.name=="Motor_vesicle_transport","cargo follows same transport rig")
            co=obj.location-Vector((0,0,2.1))
            check((co.x/1.1)**2+(co.y/.9)**2+(co.z/.9)**2<.7,"cargo centers remain within vesicle interior")

    chromatids=[o for o in scene.objects if o.type=="EMPTY" and o.name.startswith("Chromatid_")]
    if chromatids:
        scene.frame_set(scene.frame_end)
        check(len(chromatids)==8,"eight chromatids retained")
        left=[o for o in chromatids if o.matrix_world.translation.x<0]
        right=[o for o in chromatids if o.matrix_world.translation.x>0]
        check(len(left)==len(right)==4,"four daughter chromosomes on each side")
        check({o.get("chromosome_id") for o in left}=={o.get("chromosome_id") for o in right}==set(range(4)),"one chromosome identity each per daughter")
        check(bpy.data.objects["Mother_cell_cutaway"].hide_render,"mother membrane replaced after abscission")

    lead=bpy.data.objects.get("Leading_nascent_DNA_5_to_3")
    if lead:
        check(lead.data.splines[0].points[0].co.x<lead.data.splines[0].points[-1].co.x,"leading new strand points +X")
        lag=[o for o in scene.objects if o.name.startswith("Okazaki_5_to_3_")]
        check(len(lag)==4 and all(o.data.splines[0].points[0].co.x>o.data.splines[0].points[-1].co.x for o in lag),"four lagging fragments grow -X")
        scene.frame_set(1)
        initial=lead.data.bevel_factor_end
        scene.frame_set(scene.frame_end)
        check(lead.data.bevel_factor_end>initial,"leading-strand displayed length increases")

    rna=bpy.data.objects.get("Nascent_RNA_5prime_to_3prime")
    if rna:
        for frame in (1,scene.frame_end//2,scene.frame_end):
            scene.frame_set(frame)
            p=Vector(rna.data.splines[0].points[-1].co[:3])
            check((p-bpy.data.objects["RNA_polymerase_proxy"].location).length<.4,f"RNA 3-prime end stays at polymerase @ {frame}")

    peptide=bpy.data.objects.get("Nascent_peptide_N_to_C")
    if peptide:
        for frame in (1,scene.frame_end):
            scene.frame_set(frame)
            check(abs(peptide.data.splines[0].points[-1].co.z-.55)<1e-5,f"peptide C terminus stays inside ribosome proxy @ {frame}")

    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps({"file":bpy.data.filepath,"blender":bpy.app.version_string,"sample_frames":sample_frames,
                                      "checks":len(checks),"failures":failures,"result":"PASS" if not failures else "FAIL",
                                      "limits":"Prototype invariants only; visual, topology and biological review remain necessary."},indent=2)+"\n")
    print("SCENE_CHECK",len(checks),"checks",len(failures),"failures")
    if failures:
        raise RuntimeError("; ".join(failures))


if __name__=="__main__":
    main()
