#!/usr/bin/env python3
"""Сбор готовых поверхностей и диагностики для HF4-B4. Cactus не запускается.

Только стандартная библиотека Python. Исходные файлы читаются, прежние архивы
не перезаписываются. Полная серия AH3 и компактная геометрия AH1/AH2.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
import zipfile


def collect(source, destination):
    source=source.expanduser().resolve()
    if not source.is_dir():raise RuntimeError('Папка расчёта не найдена: '+str(source))
    destination=destination.expanduser().absolute()
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.is_relative_to(source):raise RuntimeError('Архив должен находиться вне папки расчёта.')
    if destination.exists():
        base=destination
        for k in range(1,10000):
            destination=base.with_name(base.stem+f'_{k:03}'+base.suffix)
            if not destination.exists():break
        else:raise RuntimeError('Не удалось выбрать свободное имя архива.')
    selected={}; vtk={}; gp={h:{} for h in [1,2,3]}
    def add(p,arc,kind):
        if p.is_symlink():raise RuntimeError('Обнаружена ссылка вместо исходного файла: '+str(p))
        if not p.is_file():return
        if p.stat().st_size>16*1024*1024 and kind=='diagnostic':return
        selected[arc]=(p,kind)
    for p in sorted(source.glob('surface03_*.vtk')):
        m=re.fullmatch(r'surface03_(\d+)\.vtk',p.name)
        if m and 9536<=int(m[1])<=11264 and int(m[1])%32==0:
            vtk[int(m[1])]=p.name;add(p,'surfaces/'+p.name,'AH3_VTK')
    if not vtk:raise RuntimeError('В выбранной папке нет ранних surface03_*.vtk. Сбор остановлен.')
    for folder in [source/'horizon_shapes',source]:
        if not folder.is_dir():continue
        for p in sorted(folder.glob('h.t*.ah*.gp')):
            m=re.fullmatch(r'h\.t(\d+)\.ah([123])\.gp',p.name)
            if not m:continue
            it,hn=int(m[1]),int(m[2]);lo=9536 if hn==3 else 0
            if lo<=it<=11264 and it%32==0:
                arc='shapes/'+p.name
                if arc in selected and selected[arc][0]!=p:
                    raise RuntimeError('Два источника одного кадра: '+p.name)
                gp[hn][it]=p.name;add(p,arc,'AH'+str(hn)+'_shape')
    for folder in [source,source/'horizon_shapes']:
        if not folder.is_dir():continue
        for pattern in ['*qlm_scalars*.asc','*qlm_state*.asc','*qlm_grid_int*.asc','*qlm_grid_real*.asc','BH_diagnostics.ah[123].gp']:
            for p in sorted(folder.glob(pattern)):
                rel=p.relative_to(source).as_posix();add(p,'diagnostics/'+rel,'diagnostic')
    for p in sorted(source.glob('*.par')):add(p,'parameters/'+p.name,'parameter')
    par=Path('/mnt/c/Users/Admin/qc0-horizon-early-geometry-t22.par')
    if par.is_file() and 'parameters/'+par.name not in selected:add(par,'parameters/'+par.name,'parameter')
    expected=set(range(9536,11265,32))
    total=sum(p.stat().st_size for p,_ in selected.values())
    if shutil.disk_usage(destination.parent).free<total+32*1024*1024:
        raise RuntimeError('Недостаточно свободного места для надёжной сборки архива.')
    print('Сбор готовых файлов. Новый расчёт не запускается.',flush=True)
    print('AH3 VTK:',len(vtk),'из 55; AH3 формы:',len(gp[3]),'из 55.',flush=True)
    print('AH1/AH2 формы:',len(gp[1]),'/',len(gp[2]),'; исходный объём:',round(total/1024**2,1),'МБ.',flush=True)
    manifest={'stage':'BH-HF4-B4-COLLECTION','created_UTC':datetime.now(timezone.utc).isoformat(),
      'source_directory':str(source),'simulation_started':False,'source_files_opened_for_write':False,
      'iterations_per_M':512,'common_birth_iteration':9536,'files':[],
      'AH3_VTK_missing_iterations':sorted(expected-set(vtk)),
      'AH3_shape_missing_iterations':sorted(expected-set(gp[3])),
      'individual_shape_missing_iterations':{str(h):sorted(set(range(0,11265,32))-set(gp[h])) for h in [1,2]},
      'frames':[{'iteration':it,'time_M':it/512,'tau_M':(it-9536)/512,
                 'vtk':vtk.get(it),'shape':gp[3].get(it)} for it in sorted(expected)]}
    with tempfile.TemporaryDirectory(prefix='hf4_b4_',dir=destination.parent) as tmp:
        temp=Path(tmp)/'collection.zip'
        with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
            for count,(arc,(p,kind)) in enumerate(sorted(selected.items()),1):
                before=p.stat();h=hashlib.sha256();size=0
                with p.open('rb') as inp,z.open(arc,'w',force_zip64=True) as output:
                    while True:
                        block=inp.read(1024*1024)
                        if not block:break
                        h.update(block);output.write(block);size+=len(block)
                after=p.stat()
                if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns) or size!=before.st_size:
                    raise RuntimeError('Исходный файл изменился во время чтения: '+str(p))
                manifest['files'].append({'archive_path':arc,'source_path':str(p),'kind':kind,'bytes':size,'sha256':h.hexdigest()})
                if count%50==0:print('Собрано',count,'из',len(selected),'файлов.',flush=True)
            manifest['common_series_complete']=not manifest['AH3_VTK_missing_iterations'] and not manifest['AH3_shape_missing_iterations']
            manifest['qlm_scalars_present']=any('qlm_scalars' in x['archive_path'] for x in manifest['files'])
            z.writestr('BH_HF4_B4_collection_manifest.json',json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
            z.writestr('READ_ME.txt','Готовые данные QC0. Новая эволюция не запускалась.\nВремена рассчитаны по подтверждённому соотношению 512 итераций на M.\nПропуски перечислены в manifest и не заполнены интерполяцией.\n')
        # Exclusive output creation preserves any file appearing meanwhile.
        with temp.open('rb') as inp,destination.open('xb') as dest:shutil.copyfileobj(inp,dest,1024*1024)
    print('Полная серия AH3:',manifest['common_series_complete'])
    print('Отдельная диагностика QLM найдена:',manifest['qlm_scalars_present'])
    print('Создан архив:',destination)
    print('Размер:',round(destination.stat().st_size/1024**2,1),'МБ')
    print('Пришлите этот архив. Исходные файлы сохранены.')
    return destination,manifest


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=Path.home()/'simulations/qc0_early_geometry_t22')
    default_dir=Path('/mnt/c/Users/Admin')
    if not default_dir.is_dir():default_dir=Path.home()
    ap.add_argument('--output',type=Path,default=default_dir/'BH_HF4_B4_collection.zip')
    args=ap.parse_args()
    try:collect(args.source,args.output)
    except Exception as e:
        print('Сбор не завершён:',e,file=sys.stderr)
        sys.exit(1)
