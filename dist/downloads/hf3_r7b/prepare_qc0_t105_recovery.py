#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BH-HF3-R7: проверенное продолжение QC0 95M→105M и сбор шва 95M.

Без флага: только подготовка. --run: подготовка, эволюция, технический сбор.
--collect: только сбор уже полученных результатов, без запуска эволюции.
Скрипт рассчитан на прежнюю WSL-среду пользователя; исходные результаты сохраняются.
"""
from __future__ import annotations
import argparse
import difflib
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
import zipfile
from pathlib import Path

VERSION = 'BH-HF3-R7-T95-T105-2026-09-15-V1'
SOURCE_SHA = '37deeb1e34a49f4fdc09dd2874a3b7f97f97d9fc0841e6c3d316d38e1f68b0d4'
START, END, START_IT, END_IT = 95., 105., 48640, 53760
MPI_PROCESSES, OMP_THREADS = 2, 2
BASENAME = 'checkpoint.chkpt'
RADII = [15, 30, 40, 50, 60]
PERIOD = 11.73680980221365
REF_KERR = {'Mf': .977461203278627, 'Jf': .648087073319352,
            'chi': .6783194994253878, 'omega': .5353401318639867,
            'alpha': .08345059699511974, 'period_M': PERIOD}
QLM_FILES = ['quasilocalmeasures-qlm_state..asc',
             'quasilocalmeasures-qlm_scalars..asc',
             'quasilocalmeasures-qlm_multipole_moments..asc']
OLD_PAR_NAME = 'qc0-horizon-multipoles-lite-t95-recover.par'
NEW_PAR_NAME = 'qc0-horizon-multipoles-lite-t105-recover.par'
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n',encoding='utf-8')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def source_text(path):
    # Смена Windows CRLF на LF не меняет параметры и допускается.
    text = path.read_text(encoding='utf-8')
    actual = hashlib.sha256(text.encode('utf-8')).hexdigest()
    if actual != SOURCE_SHA:
        raise RuntimeError('Исходный файл t95 отличается от подтверждённого R6A. '
                           'Новый расчёт не запущен; отличие включено в диагностический архив.')
    return text


def get_parameter(text, key):
    pattern = r'^\s*'+re.escape(key)+r'\s*=\s*("[^"\n]*"|[^#\n]+)'
    hits = re.findall(pattern,text,re.MULTILINE|re.IGNORECASE)
    if len(hits) != 1:
        raise RuntimeError(f'{key}: ожидалась ровно одна запись, найдено {len(hits)}')
    return hits[0].strip().strip('"')


def build_parameter(text, old_dir, new_dir):
    for key, expected in [('Cactus::cctk_final_time','95'),('Cactus::terminate','time'),
                          ('IO::recover','auto'),('IO::recover_file',BASENAME),
                          ('IO::recover_and_remove','no'),
                          ('IO::truncate_files_after_recovering','no'),
                          ('IO::checkpoint_on_terminate','yes')]:
        if get_parameter(text,key).lower() != expected.lower():
            raise RuntimeError(f'Неожиданное значение {key}')
    changes = {'Cactus::cctk_final_time':'105',
               'IO::out_dir': f'"{new_dir}"',
               'IO::checkpoint_dir': f'"{new_dir}"',
               'IO::recover_dir': f'"{old_dir}"'}
    result = text
    for key, value in changes.items():
        pattern=r'^(\s*'+re.escape(key)+r'\s*=\s*)("[^"\n]*"|[^#\n]+)'
        result,n=re.subn(pattern,lambda m:m.group(1)+value,result,
                         flags=re.MULTILINE|re.IGNORECASE)
        if n != 1:
            raise RuntimeError(f'Неоднозначная замена {key}')
    # Обратная замена доказывает неизменность всех остальных байтов тела .par.
    restored = result
    for key in changes:
        pattern=r'^(\s*'+re.escape(key)+r'\s*=\s*)("[^"\n]*"|[^#\n]+)'
        original = re.search(pattern,text,re.MULTILINE|re.IGNORECASE).group(2)
        restored=re.sub(pattern,lambda m:m.group(1)+original,restored,
                        flags=re.MULTILINE|re.IGNORECASE)
    if restored != text:
        raise RuntimeError('Обнаружено изменение за пределами четырёх разрешённых параметров')
    header=(f'# {VERSION}\n# CURRENT RECOVERY: t=95M, it=48640 -> t=105M, it=53760\n'
            '# All older generated comments below record previous stages.\n'
            '# Физические параметры и параметры сетки сохранены из подтверждённого t95.\n\n')
    return header+result, changes


def checkpoint_inventory(directory, iteration, require_latest=False):
    pattern=re.compile(re.escape(BASENAME)+r'\.it_(\d+)\.file_(\d+)\.h5$')
    found=[]
    all_iterations=[]
    for p in directory.glob(BASENAME+'.it_*.file_*.h5'):
        m=pattern.fullmatch(p.name)
        if m:
            it,part=map(int,m.groups());all_iterations.append(it)
            if it == iteration:
                found.append((part,p))
    if [part for part,_ in sorted(found)] != [0,1]:
        raise RuntimeError(f'Для it={iteration} нужны ровно части file_0.h5 и file_1.h5: {directory}')
    if require_latest and max(all_iterations) != iteration:
        raise RuntimeError('В каталоге есть более поздняя контрольная точка. '
                           'Запуск из 95M остановлен, чтобы не повторять уже выполненное.')
    records=[]
    for part,p in sorted(found):
        st=p.stat()
        if st.st_size < 512:
            raise RuntimeError(f'Пустая или неполная контрольная точка: {p}')
        with p.open('rb') as f:
            magic=f.read(8)
        if magic != b'\x89HDF\r\n\x1a\n':
            raise RuntimeError(f'Не распознан заголовок HDF5: {p}')
        records.append({'path':str(p),'part':part,'bytes':st.st_size,
                        'mtime_ns':st.st_mtime_ns,'inode':st.st_ino,
                        'header_is_hdf5':True,
                        'check_scope':'Имя, размер и заголовок. Полное содержимое проверяет загрузчик Cactus.'})
    return records


def load_array(path):
    import numpy as np
    a=np.loadtxt(path,comments='#',ndmin=2)
    if a.size == 0:
        raise RuntimeError(f'Пустой ряд: {path}')
    return a


def column_map(path):
    with path.open(encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('# data columns:'):
                return {name:int(j)-1 for j,name in re.findall(r'(\d+):([^\s]+)',line)}
    return {}


def require_old_end(directory):
    import numpy as np
    for r in RADII:
        p=directory/f'mp_Psi4_l2_m2_r{r:.2f}.asc'
        a=load_array(p)
        if a.shape[1] != 3 or not np.isclose(a[-1,0],START,atol=1e-10,rtol=0):
            raise RuntimeError(f'Ряд {p.name} не заканчивается на 95M')
        if not np.all(np.isfinite(a[-1])):
            raise RuntimeError(f'Некорректные числа в конце {p.name}')
    state=directory/QLM_FILES[0]; scalars=directory/QLM_FILES[1]
    a,b=load_array(state),load_array(scalars)
    ca,cb=column_map(state),column_map(scalars)
    if a[-1,ca['qlm_have_valid_data[2]']] <= .5 or abs(b[-1,cb['qlm_time[2]']]-START)>1e-10:
        raise RuntimeError('Последняя запись общего горизонта slot 2 не подтверждает 95M')


def paths(args):
    root=args.windows_root.resolve()
    return {'root':root, 'source':args.source_par or root/OLD_PAR_NAME,
            'old':args.recovery_dir.resolve(), 'new':args.new_output_dir.resolve(),
            'cactus':args.cactus_dir.resolve(), 'par':root/NEW_PAR_NAME,
            'manifest':root/'qc0-horizon-multipoles-lite-t105-recovery_manifest.json',
            'diff':root/'BH_HF3_R7_parameter_changes.diff',
            'log':root/'ET_QC0_lite_t105_recover.log',
            'status':root/'BH_HF3_R7_run_status.json',
            'summary':root/'BH_HF3_R7A_summary.json',
            'bundle':root/'BH_HF3_R7A_output.zip',
            'diagnostic':root/'BH_HF3_R7_diagnostic.zip'}


def write_same_or_new(path, text):
    if path.exists() and path.read_text(encoding='utf-8') != text:
        raise RuntimeError(f'Файл уже существует с другим содержимым: {path}')
    if not path.exists():
        with path.open('x',encoding='utf-8',newline='\n') as f:
            f.write(text)


def prepare(p):
    if not p['root'].is_dir():
        raise RuntimeError(f'Не найдена папка Windows: {p["root"]}. Запустите файл в прежнем WSL.')
    text=source_text(p['source'])
    if not p['old'].is_dir():
        raise RuntimeError(f'Не найден каталог контрольной точки 95M: {p["old"]}')
    if p['new'] == p['old'] or p['new'] in p['old'].parents or p['old'] in p['new'].parents:
        raise RuntimeError('Каталоги старого и нового расчёта должны быть отдельными')
    if p['new'].exists() and any(p['new'].iterdir()):
        raise RuntimeError('Каталог t105 уже содержит результаты. Для их сбора используйте --collect; '
                           'повторная эволюция не запущена.')
    if p['log'].exists():
        raise RuntimeError('Журнал t105 уже существует. Используйте --collect для проверки имеющегося запуска.')
    checkpoints=checkpoint_inventory(p['old'],START_IT,require_latest=True)
    require_old_end(p['old'])
    executable=p['cactus']/'exe/cactus_sim'
    if not executable.is_file() or not os.access(executable,os.X_OK):
        raise RuntimeError(f'Не найдена прежняя исполняемая программа: {executable}')
    if not shutil.which('mpirun'):
        raise RuntimeError('В прежней WSL-среде не найден mpirun')
    if not Path('/usr/bin/time').is_file():
        raise RuntimeError('Не найден /usr/bin/time для записи результата запуска')
    if shutil.which('pgrep'):
        active=subprocess.run(['pgrep','-u',str(os.getuid()),'-x','cactus_sim'],
                              capture_output=True,text=True)
        if active.returncode == 0:
            raise RuntimeError('У этого пользователя уже работает cactus_sim; второй расчёт не запущен.')
    anchor=p['new'].parent
    while not anchor.exists():
        anchor=anchor.parent
    free=shutil.disk_usage(anchor).free
    minimum=int(sum(c['bytes'] for c in checkpoints)*1.2+512*1024**2)
    if free < minimum:
        raise RuntimeError(f'Для новых контрольных точек требуется около {minimum/1024**3:.2f} GiB; '
                           f'доступно {free/1024**3:.2f} GiB. Старые файлы не изменены.')
    newtext,changes=build_parameter(text,p['old'],p['new'])
    write_same_or_new(p['par'],newtext)
    diff=''.join(difflib.unified_diff(text.splitlines(True),newtext.splitlines(True),
                                    fromfile=OLD_PAR_NAME,tofile=NEW_PAR_NAME))
    write_same_or_new(p['diff'],diff)
    manifest={'stage':'BH-HF3-R7','version':VERSION,'state':'PREPARED_NOT_RUN',
              'source_parameter_file':str(p['source']),'source_normalized_sha256':SOURCE_SHA,
              'generated_parameter_file':str(p['par']),'generated_sha256':sha(p['par']),
              'changed_parameters':changes,'all_other_parameter_text_preserved':True,
              'recovery_directory':str(p['old']),'new_output_directory':str(p['new']),
              'recovery_iteration':START_IT,'expected_start_time_M':START,
              'target_final_time_M':END,'expected_final_iteration_if_step_unchanged':END_IT,
              'checkpoint_files':checkpoints,'checkpoint_content_fully_read_by_preparer':False,
              'recovery_mode':'auto','latest_checkpoint_iteration_verified':START_IT,
              'recover_and_remove':False,'truncate_files_after_recovering':False,
              'working_directory':str(p['cactus']),'executable':str(executable),
              'mpi_processes':MPI_PROCESSES,'omp_threads_per_process':OMP_THREADS,
              'launch_evidence':'Журнал подтверждённого R6: mpirun -n 2 ./exe/cactus_sim; 2 threads per process.',
              'free_space_bytes':free,'minimum_space_precheck_bytes':minimum,
              'QC0_numerical_Kerr_reference':REF_KERR,
              'scientific_goal':{'r40_third_period_end_M':69+3*PERIOD,
                                 'r50_second_period_end_M':79.25+2*PERIOD,
                                 'r60_first_period_end_M':89.75+PERIOD},
              'program_sha256':sha(Path(__file__).resolve())}
    dump(p['manifest'],manifest)
    print('Подготовка завершена: 95M / it=48640 → 105M / it=53760.',flush=True)
    print('Все параметры, кроме конечного времени и трёх путей, сохранены.',flush=True)
    print('Проверены обе части контрольной точки и окончание прежних рядов на 95M.',flush=True)
    print('Файл параметров:',p['par'],flush=True)
    return manifest


def monitor_line(line, state):
    recovery_file=re.search(r'Recovering parameters from checkpoint file.*\.it_(\d+)\.file_',line)
    if recovery_file and int(recovery_file.group(1)) != START_IT:
        raise RuntimeError('Cactus выбрал контрольную точку другой итерации')
    restart=re.search(r'restarting simulation .* at iteration (\d+) \(simulation time ('+NUMBER+r')\)',line)
    if restart:
        it,t=int(restart.group(1)),float(restart.group(2))
        if it != START_IT or abs(t-START)>1e-9:
            raise RuntimeError(f'Восстановление началось с it={it}, t={t}; ожидалось it=48640, t=95')
        state['restart_confirmations']+=1
    progress=re.match(r'^\s*(\d+)\s+('+NUMBER+r')\s*\|',line)
    if progress:
        it,t=int(progress.group(1)),float(progress.group(2))
        if not state['restart_confirmations']:
            raise RuntimeError('Пошла эволюция без подтверждения восстановления из 95M')
        if it < START_IT or it > END_IT or t < START-1e-9 or t > END+1e-9:
            raise RuntimeError(f'Эволюция вышла за ожидаемые пределы: it={it}, t={t}')
        if state['last_iteration'] is not None and it < state['last_iteration']:
            raise RuntimeError('Номер итерации уменьшился')
        state['last_iteration'],state['last_time_M']=it,t
        if state['last_printed_integer_M'] != int(t):
            state['last_printed_integer_M']=int(t)
            return f'Идёт расчёт: t={t:.3f}M, it={it}'
    return None


def stop_own_process(proc):
    if proc.poll() is not None:
        return
    for sig,timeout in [(signal.SIGINT,12),(signal.SIGTERM,8),(signal.SIGKILL,5)]:
        try:
            os.killpg(proc.pid,sig)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=timeout)
            return
        except subprocess.TimeoutExpired:
            pass


def run_evolution(p, manifest):
    # Между подготовкой и запуском повторно проверяются только критические входы.
    if sha(p['par']) != manifest['generated_sha256']:
        raise RuntimeError('Подготовленный .par изменился перед запуском')
    before=checkpoint_inventory(p['old'],START_IT,require_latest=True)
    if before != manifest['checkpoint_files']:
        raise RuntimeError('Контрольная точка изменилась между подготовкой и запуском')
    p['new'].mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['OMP_NUM_THREADS']=str(OMP_THREADS)
    command=['/usr/bin/time','-v',shutil.which('mpirun'),'-n',str(MPI_PROCESSES),
             './exe/cactus_sim',str(p['par'])]
    state={'stage':'BH-HF3-R7','version':VERSION,'status':'RUNNING',
           'command':command,'working_directory':str(p['cactus']),
           'OMP_NUM_THREADS':OMP_THREADS,'start_unix_seconds':time.time(),
           'restart_confirmations':0,'last_iteration':None,'last_time_M':None,
           'last_printed_integer_M':None,'process_exit_code':None,'guard_error':None}
    dump(p['status'],state)
    print('Запускается продолжение. Полный журнал:',p['log'],flush=True)
    print('Оставьте это окно WSL открытым до сообщения о завершении.',flush=True)
    started=time.monotonic();proc=None
    with p['log'].open('x',encoding='utf-8') as log:
        try:
            proc=subprocess.Popen(command,cwd=p['cactus'],env=env,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                  text=True,encoding='utf-8',errors='replace',bufsize=1,
                                  start_new_session=True)
            for line in proc.stdout:
                log.write(line);log.flush()
                msg=monitor_line(line,state)
                if msg:
                    print(msg,flush=True);dump(p['status'],state)
            state['process_exit_code']=proc.wait()
        except BaseException as e:
            state['guard_error']=str(e) or type(e).__name__
            if proc is not None:
                stop_own_process(proc)
                state['process_exit_code']=proc.poll()
            raise
        finally:
            state['elapsed_seconds']=time.monotonic()-started
            state['status']='EXITED' if state['process_exit_code']==0 and not state['guard_error'] else 'STOPPED_OR_FAILED'
            dump(p['status'],state)
    if state['process_exit_code'] != 0:
        raise RuntimeError(f'Cactus завершился с кодом {state["process_exit_code"]}; журнал сохранён')
    if checkpoint_inventory(p['old'],START_IT,require_latest=True) != before:
        raise RuntimeError('После запуска изменились исходные контрольные точки 95M')
    print('Эволюция завершилась. Проверяю шов 95M и собираю данные.',flush=True)


def validate_grid(times,start,end,step):
    import numpy as np
    expected=np.arange(round((end-start)/step)+1)*step+start
    return bool(times.shape==expected.shape and np.allclose(times,expected,atol=1e-10,rtol=0))


def audit_mode(old,new):
    import numpy as np
    a,b=load_array(old),load_array(new)
    oldseam=np.flatnonzero(np.isclose(a[:,0],START,atol=1e-10,rtol=0))
    newseam=np.flatnonzero(np.isclose(b[:,0],START,atol=1e-10,rtol=0))
    shape= a.shape[1]==b.shape[1]==3
    finite=bool(np.all(np.isfinite(b)))
    exact=bool(shape and len(oldseam)==1 and len(newseam)==1 and
               np.array_equal(a[oldseam[0]],b[newseam[0]]))
    grid=shape and validate_grid(b[:,0],START,END,.25)
    return {'file':new.name,'samples':len(b),'finite':finite,
            'complete_grid_95_105':grid,'seam_exact':exact,'pass':bool(finite and grid and exact)}


def audit_qlm(old,new):
    import numpy as np
    a,b=load_array(old),load_array(new)
    ma,mb=column_map(old),column_map(new)
    fields=sorted(k for k in ma if k.endswith('[2]'))
    complete=bool(fields and all(k in mb for k in fields))
    grid=validate_grid(b[:,8],START,END,.0625)
    ia=np.flatnonzero(np.isclose(a[:,8],START,atol=1e-10,rtol=0))
    ib=np.flatnonzero(np.isclose(b[:,8],START,atol=1e-10,rtol=0))
    exact=False;finite=False;maxdiff=None
    if complete:
        x,y=a[:,[ma[k] for k in fields]],b[:,[mb[k] for k in fields]]
        finite=bool(np.all(np.isfinite(y)))
        if len(ia)==1 and len(ib)==1:
            exact=bool(np.array_equal(x[ia[0]],y[ib[0]]))
            if np.all(np.isfinite(x[ia[0]])) and np.all(np.isfinite(y[ib[0]])):
                maxdiff=float(np.max(abs(x[ia[0]]-y[ib[0]])))
    valid=True
    if 'qlm_have_valid_data[2]' in mb:
        valid=bool(np.all(b[:,mb['qlm_have_valid_data[2]']]>0.5))
    if 'qlm_time[2]' in mb:
        valid=valid and bool(np.allclose(b[:,mb['qlm_time[2]']],b[:,8],atol=1e-10,rtol=0))
    return {'file':new.name,'slot':2,'fields':fields,'samples':len(b),
            'headers_compatible':complete,'slot_fields_finite':finite,
            'complete_grid_95_105':grid,'slot_valid_or_current':valid,
            'seam_exact':exact,'maximum_seam_difference':maxdiff,
            'pass':bool(complete and finite and grid and valid and exact)}


def collect(p):
    if not p['new'].is_dir():
        raise RuntimeError('Результатов t105 пока нет')
    records=[];errors=[]
    for r in RADII:
        for ell in range(5):
            for m in range(-ell,ell+1):
                name=f'mp_Psi4_l{ell}_m{m}_r{r:.2f}.asc'
                try:
                    records.append(dict(audit_mode(p['old']/name,p['new']/name),radius_M=r))
                except Exception as e:
                    errors.append({'file':name,'error':str(e)})
    qlm=[]
    for name in QLM_FILES:
        try:
            qlm.append(audit_qlm(p['old']/name,p['new']/name))
        except Exception as e:
            errors.append({'file':name,'error':str(e)})
    log=p['log'].read_text(encoding='utf-8',errors='replace') if p['log'].exists() else ''
    status=json.loads(p['status'].read_text()) if p['status'].exists() else {}
    log_checks={'restarted_it48640_t95':bool(re.search(r'restarting simulation .* at iteration 48640 \(simulation time 95(?:\.0*)?\)',log)),
                'termination_time_105':'Terminating due to cctk_final_time at t = 105.000000' in log,
                'termination_checkpoint_53760':'Dumping termination checkpoint at iteration 53760, simulation time 105' in log,
                'native_exit_0':bool(re.search(r'^\s*Exit status:\s*0\s*$',log,re.MULTILINE)),
                'wrapper_exit_0':status.get('process_exit_code')==0 and not status.get('guard_error')}
    checkpoints=[]
    try:
        checkpoints=checkpoint_inventory(p['new'],END_IT)
    except Exception as e:
        errors.append({'checkpoint':str(e)})
    old_unchanged=False
    try:
        manifest=json.loads(p['manifest'].read_text())
        old_unchanged=checkpoint_inventory(p['old'],START_IT,require_latest=True)==manifest['checkpoint_files']
    except Exception as e:
        errors.append({'source_checkpoint':str(e)})
    passed=bool(not errors and len(records)==125 and all(r['pass'] for r in records) and
                len(qlm)==3 and all(r['pass'] for r in qlm) and checkpoints and
                all(log_checks.values()) and old_unchanged)
    summary={'stage':'BH-HF3-R7A','version':VERSION,'technical_pass':passed,
             'classification':'BH_HF3_R7A_T105_SEAM95_DATA_READY' if passed else 'BH_HF3_R7A_REVIEW_REQUIRED',
             'run_log_checks':log_checks,'new_checkpoint_files':checkpoints,
             'source_checkpoints_unchanged':old_unchanged,'Psi4_modes':records,'qlm_slot2':qlm,
             'errors':errors,'scientific_fits_performed':False,
             'scientific_goal_end_M':{'r40_third_period':69+3*PERIOD,
                                      'r50_second_period':79.25+2*PERIOD,
                                      'r60_first_period':89.75+PERIOD},
             'QC0_numerical_Kerr_reference_from_R6B':REF_KERR,
             'next_step':'Проверка поздних периодов по данным t105, с одинаковым комплексным критерием для фиксированной и свободной моделей.'}
    dump(p['summary'],summary)
    report=p['root']/'BH_HF3_R7A_report.txt'
    report.write_text('BH-HF3-R7A\n'+summary['classification']+'\n'
                      f'Технический проход: {passed}\n'
                      f'Контроль мод Psi4: {sum(r["pass"] for r in records)}/125\n'
                      f'Контроль общего горизонта: {sum(r["pass"] for r in qlm)}/3\n'
                      f'Исходные контрольные точки сохранены: {old_unchanged}\n'
                      'Физическое сравнение кольцевания выполняется после чтения этого пакета.\n',encoding='utf-8')
    # Сначала временный новый ZIP; замена только нашего предыдущего технического пакета.
    temporary=p['bundle'].with_suffix('.zip.tmp')
    with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED) as z:
        for key in ['source','par','manifest','diff','log','status','summary']:
            if p[key].is_file():
                z.write(p[key],'provenance/'+p[key].name)
        z.write(report,report.name)
        z.write(Path(__file__).resolve(),'provenance/'+Path(__file__).name)
        for directory,label in [(p['old'],'t70_t95'),(p['new'],'t95_t105')]:
            selected=list(directory.glob('mp_Psi4_l*_m*_r*.asc'))
            selected+=list(directory.glob('quasilocalmeasures-*.asc'))
            selected+=list(directory.glob('sphericalsurface-*.asc'))
            for f in sorted(set(selected)):
                if f.is_file():
                    z.write(f,label+'/'+f.name)
    os.replace(temporary,p['bundle'])
    print(report.read_text(),flush=True)
    print('Пришлите один архив:',p['bundle'],flush=True)
    return summary


def diagnostic(p, error):
    if not p['root'].is_dir():
        return None
    info=p['root']/'BH_HF3_R7_diagnostic.json'
    dump(info,{'stage':'BH-HF3-R7','version':VERSION,'error':str(error),
               'run_requested':('--run' in sys.argv),'traceback':traceback.format_exc(),
               'original_results_deleted':False})
    with zipfile.ZipFile(p['diagnostic'],'w',zipfile.ZIP_DEFLATED) as z:
        z.write(info,info.name)
        for key in ['source','par','manifest','diff','status','summary']:
            if p[key].is_file():
                z.write(p[key],p[key].name)
        if p['log'].is_file():
            # Журнал этого запуска полезен целиком; контрольные точки в ZIP не входят.
            z.write(p['log'],p['log'].name)
    return p['diagnostic']


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group()
    mode.add_argument('--run',action='store_true')
    mode.add_argument('--collect',action='store_true')
    ap.add_argument('--windows-root',type=Path,default=Path('/mnt/c/Users/Admin'))
    ap.add_argument('--source-par',type=Path)
    ap.add_argument('--recovery-dir',type=Path,default=Path('/home/alex/simulations/qc0_lite_t95_recover'))
    ap.add_argument('--new-output-dir',type=Path,default=Path('/home/alex/simulations/qc0_lite_t105_recover'))
    ap.add_argument('--cactus-dir',type=Path,default=Path('/home/alex/EinsteinToolkit_2026_05/Cactus'))
    args=ap.parse_args();p=paths(args)
    lock=None
    try:
        import numpy  # Уже используется в прежней среде gw_env; проверяется до запуска.
        import fcntl
        if not p['root'].is_dir():
            raise RuntimeError('Нужна прежняя WSL-среда с /mnt/c/Users/Admin')
        lock=(p['root']/'BH_HF3_R7.lock').open('a+')
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print('Этот сценарий уже работает; второй экземпляр не запущен.',flush=True)
            return 2
        if args.collect:
            result=collect(p)
            return 0 if result['technical_pass'] else 2
        manifest=prepare(p)
        if args.run:
            run_evolution(p,manifest)
            result=collect(p)
            return 0 if result['technical_pass'] else 2
        print('Для выполнения подготовленного продолжения добавьте --run.',flush=True)
        return 0
    except KeyboardInterrupt as e:
        print('Остановка по Ctrl+C. Точка 95M и прежние результаты сохранены.',flush=True)
        archive=diagnostic(p,e)
        if archive: print('Диагностический архив:',archive,flush=True)
        return 130
    except Exception as e:
        print('Остановлено:',e,flush=True)
        archive=diagnostic(p,e)
        if archive: print('Пришлите диагностический архив:',archive,flush=True)
        return 2
    finally:
        if lock is not None:
            lock.close()


if __name__=='__main__':
    sys.exit(main())
