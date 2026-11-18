from pathlib import Path
import json,sys,shutil,hashlib,csv
import fitz
def build(historical,current,dest):
 dest.mkdir(exist_ok=True,parents=True)
 def match(rel):
  return historical/rel,current/rel
 def add(book,toc,rel,title,nums=None):
  left,right=match(rel)
  with fitz.open(left) as a,fitz.open(right) as b:
   assert len(a)==len(b),(rel,len(a),len(b))
   numbers=list(range(len(a))) if nums is None else nums
   toc.append([1,title,len(book)+1])
   for k in numbers:
    page=book.new_page(width=1200,height=500);page.insert_text((22,23),title+'   |   '+str(k+1)+' / '+str(len(a)),fontsize=14)
    for t,source,xpos in [('Published Run 25521',a,18),('Current 222-cluster indexing',b,608)]:
     page.insert_text((xpos,43),t,fontsize=11,color=(.10,.31,.43));page.show_pdf_page(fitz.Rect(xpos,51,xpos+568,477),source,k)
    page.insert_text((20,492),'Unweighted Run 25521. Each panel retains its native axes, scale and fitted population. Historical: 73,188; current: 73,989. Available-event delta panels use 106,058.',fontsize=7)
 rels=[('plots/25521_available_delta_histograms.pdf','Available delta distributions'),('plots/25521_delta_histograms.pdf','Fitted delta distributions'),('unweighted/plots/comparison/25521_delta_vs_ytar_before_after.pdf','Starting and fitted matrices: delta versus ytar')]
 for side in ['old','new']:
  folder=current/f'unweighted/plots/{side}'
  rels.append((f'unweighted/plots/fit_opt_matrix_25521_unweighted_{side}_matrix_diff.pdf',f'{side.title()} matrix: automatic summary'))
  for path in sorted(folder.glob('*.pdf')):rels.append((str(path.relative_to(current)),side+' / '+path.stem))
  for path in sorted((folder/'nps_legacy').glob('*.pdf')):
   if not path.name.startswith('foil_'):rels.append((str(path.relative_to(current)),side+' / NPS / '+path.stem))
  for path in sorted(folder.glob('foil_*/all_delta_slices.pdf')):rels.append((str(path.relative_to(current)),side+' / '+path.parent.name+' / six delta slices'))
 book=fitz.open();toc=[]
 for rel,name in rels:add(book,toc,rel,name)
 book.set_toc(toc);book.save(dest/'unweighted_complete_atlas.pdf',garbage=4,deflate=True);book.close()
 quick=fitz.open();toc=[]
 for rel,name in rels[:3]:add(quick,toc,rel,name)
 for f,foil in enumerate(['foil_0_z_plus10','foil_1_z_0','foil_2_z_minus10']):
  add(quick,toc,f'unweighted/plots/new/{foil}/delta_04_plus6_to_plus10.pdf',foil+' / 6 <= delta < 10 / sieve and three residuals',[4])
  add(quick,toc,'unweighted/plots/new/nps_legacy/25521_new_xpdiff.pdf',foil+' / 6 <= delta < 10 / xptar residual by hole',[f*6+4])
 quick.set_toc(toc);quick.save(dest/'unweighted_highlights.pdf',garbage=4,deflate=True)
 print('atlases',len(quick),len(rels))
if __name__=='__main__':build(*(Path(s).resolve() for s in sys.argv[1:]))
