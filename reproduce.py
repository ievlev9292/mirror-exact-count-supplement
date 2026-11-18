"""Reproduce the manuscript figures and compile the theorem proof.

Run `python reproduce.py`; PDFs are written to figures/ and proof/. Temporary
plot tables and TeX/BibTeX files stay in a separate working directory. The
parameters below match the accompanying manuscript; no data download is used.
"""
import os
for variable in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(variable, '1')
from pathlib import Path
import argparse
import math
import shutil
import subprocess
import tempfile
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT/'figures'
PEAK_SPEED = 0.97
NEGATIVE_SPEED = 0.65
CURVES = ('m1', 'm2', 'm3', 'bounce', 'asymmetric', 'elliptic')


def bisect(function, lower, upper):
    """Solve a bracketed increasing equation at the current mp precision."""
    lower, upper = mp.mpf(lower), mp.mpf(upper)
    if function(lower) > 0 or function(upper) < 0:
        raise ValueError('The increasing equation is not bracketed.')
    for _ in range(4*mp.mp.dps+20):
        middle = (lower+upper)/2
        if middle == lower or middle == upper:
            break
        if function(middle) < 0:
            lower = middle
        else:
            upper = middle
    return (lower+upper)/2


def elliptic_parameters(speed):
    """Calibrate L/tau from the actual maximum speed, not its large-L limit.

    At the positive speed extremum h=z*/L lies in (0,1/4). Differentiating
    t^2=z^2+tau^2*z/(L-z) gives
    4*h^3=(1-V^2)*(1-3*h)^2, followed by the L/tau formula below.
    """
    speed = mp.mpf(speed)
    if not 0 <= speed < 1:
        raise ValueError('The peak speed must satisfy 0<=V<1.')
    if speed == 0:
        return mp.mpf(0), mp.mpf(0), mp.mpf(0)
    peak_fraction = bisect(lambda h:4*h**3-(1-speed**2)*(1-3*h)**2, 0, mp.mpf(1)/4)
    excursion = mp.sqrt((1-4*peak_fraction)/(4*peak_fraction**3*(1-peak_fraction)))
    peak_time = (1-2*peak_fraction)/(2*mp.sqrt(peak_fraction*(1-peak_fraction)))
    return excursion, peak_time, peak_fraction


def count_coefficient(order):
    """Positive cotangent coefficient C_m, for an integer m>=1."""
    if isinstance(order, bool) or int(order) != order or order < 1:
        raise ValueError('The peak-order index must be a positive integer.')
    order = int(order)
    return mp.fsum(mp.cot(mp.pi*p/(2*order+1)) for p in range(1,order+1))/(mp.pi*order)


def total_energy(order, speed):
    """Exact two-sided tau*E for a polynomial flight, 0<=V<1."""
    speed = mp.mpf(speed)
    if not 0 <= speed < 1:
        raise ValueError('The maximum speed must satisfy 0<=V<1.')
    if speed == 0:
        return mp.mpf(0)
    return ((1-mp.mpf(1)/(2*order))/(48*mp.sin(mp.pi/(2*order)))
            *((1-speed)**(-1-mp.mpf(1)/(2*order))
              +(1+speed)**(-1-mp.mpf(1)/(2*order))
              +(2*order/speed)*((1+speed)**(-mp.mpf(1)/(2*order))
                               -(1-speed)**(-mp.mpf(1)/(2*order)))))


def worldline(name, coordinate, excursion=0.0):
    """Return t/tau,z/tau and derivatives with respect to the sample coordinate.

    The coordinate is only an integration/sampling variable. All curves use
    their defining origins. For the elliptic trip z=L*x^2/(1+x^2), while
    t/tau=x*sqrt(1+(L/tau)^2*x^2/(1+x^2)^2). Thus it turns at t=z=0
    and approaches z=L at both temporal ends.
    """
    if name in ('m1','m2','m3'):
        order = int(name[1:])
        return (coordinate+coordinate**(2*order+1)/(2*order+1),
                -PEAK_SPEED*coordinate, 1+coordinate**(2*order), -PEAK_SPEED)
    if name in ('bounce','asymmetric'):
        other_speed = PEAK_SPEED if name == 'bounce' else NEGATIVE_SPEED
        position = ((PEAK_SPEED-other_speed)*coordinate
                    +math.sqrt(PEAK_SPEED*other_speed)*coordinate**2)
        position_derivative = PEAK_SPEED-other_speed+2*math.sqrt(PEAK_SPEED*other_speed)*coordinate
        return coordinate+coordinate**3/3, position, 1+coordinate**2, position_derivative
    if name == 'elliptic':
        square = 1+coordinate**2
        factor = math.sqrt(1+(excursion*coordinate/square)**2)
        return (coordinate*factor, excursion*coordinate**2/square,
                (1+2*excursion**2*coordinate**2/square**3)/factor,
                2*excursion*coordinate/square**2)
    raise ValueError('Unknown curve: '+name)


def compactify(time, position):
    """X=(atan(t+z)-atan(t-z))/pi; T=(atan(t+z)+atan(t-z))/pi."""
    advanced = math.atan(time+position)
    retarded = math.atan(time-position)
    return (advanced-retarded)/math.pi, (advanced+retarded)/math.pi


def check_results():
    """Small normalization and geometry controls; not a proof of the theorem."""
    with mp.workdps(60):
        assert abs(bisect(lambda x:x-mp.mpf(3)/16,0,1)-mp.mpf(3)/16) < mp.mpf('1e-55')
        try:
            bisect(lambda x:x-2,0,1)
        except ValueError:
            pass
        else:
            raise AssertionError('An invalid bracket was accepted.')
        assert worldline('elliptic',2,0) == (2.0,0.0,1.0,0.0)
        planted_excursion,planted_time,planted_fraction = elliptic_parameters(7/mp.sqrt(50))
        assert abs(planted_fraction-mp.mpf(1)/8) < mp.mpf('1e-55')
        assert abs(planted_excursion-16*mp.sqrt(mp.mpf(2)/7)) < mp.mpf('1e-54')
        assert abs(planted_time-3/mp.sqrt(7)) < mp.mpf('1e-55')
        assert abs(count_coefficient(1)-1/(mp.pi*mp.sqrt(3))) < mp.mpf('1e-55')
        assert abs(count_coefficient(2)-mp.sqrt((5+mp.sqrt(5))/10)/mp.pi) < mp.mpf('1e-55')
        for value in ('0.2','0.7','0.97'):
            speed = mp.mpf(value)
            gamma = 1/mp.sqrt(1-speed**2)
            cubic = (gamma-1)*(2*gamma+3)*gamma**mp.mpf('1.5')/(48*mp.sqrt(2*(gamma+1)))
            assert abs(total_energy(1,speed)-cubic) < mp.mpf('1e-52')
        excursion,peak_time,fraction = elliptic_parameters(str(PEAK_SPEED))
        peak_coordinate = float(mp.sqrt(fraction/(1-fraction)))
        _,_,time_derivative,position_derivative = worldline('elliptic',peak_coordinate,float(excursion))
        assert abs(position_derivative/time_derivative-PEAK_SPEED) < 2e-14
    assert compactify(0,0) == (0,0)
    assert abs(compactify(1,0)[1]-0.5) < 1e-14
    for name in CURVES:
        for coordinate in (-1000,-2,-1,0,1,2,1000):
            time,position,time_derivative,position_derivative = worldline(name,coordinate,float(excursion))
            assert time_derivative > abs(position_derivative)
            x_value,t_value = compactify(time,position)
            assert abs(x_value)+abs(t_value) <= 1+1e-14
    return float(excursion),float(peak_time)


def write_table(folder,name,columns):
    keys = list(columns)
    if len({len(column) for column in columns.values()}) != 1:
        raise ValueError('Unequal column lengths in '+name)
    text = ' '.join(keys)+'\n'+'\n'.join(
        ' '.join(format(value,'.16g') for value in row)
        for row in zip(*columns.values()))+'\n'
    (folder/name).write_text(text,encoding='utf-8',newline='')


def make_data(work):
    data = work/'data'
    data.mkdir(parents=True,exist_ok=True)
    excursion,peak_time = check_results()
    with mp.workdps(60):
        coefficients = np.array([float(count_coefficient(order)) for order in range(1,8)])
        write_table(data,'coefficients.dat',{'m':list(range(1,8)),'ratio':coefficients/coefficients[0]})
        gamma_values = np.geomspace(1.02,10000,180)
        energies = {order:[] for order in (2,3,4)}
        for gamma_value in gamma_values:
            speed = mp.sqrt(1-1/mp.mpf(str(gamma_value))**2)
            cubic_energy = total_energy(1,speed)
            for order in energies:
                energies[order].append(float(total_energy(order,speed)/cubic_energy))
        write_table(data,'energies.dat',{'gamma':gamma_values,
                    **{'m'+str(order):values for order,values in energies.items()}})
        _,_,fraction = elliptic_parameters(str(PEAK_SPEED))
        peak_coordinate = float(mp.sqrt(fraction/(1-fraction)))
    coordinate_values = np.linspace(-4,4,5001)
    for name,table_name in zip(CURVES,('flight1.dat','flight2.dat','flight3.dat','return0.dat','return1.dat','elliptic.dat')):
        samples = [coordinate for index,coordinate in enumerate(coordinate_values) if index%3==0]
        if name == 'elliptic':
            samples = sorted(set(np.linspace(-4,4,10001)[::3].tolist()
                                 +[-peak_coordinate,0.0,peak_coordinate]))
        rows = []
        for coordinate in samples:
            time,_,time_derivative,position_derivative = worldline(name,coordinate,excursion)
            if -3.6 < time < 3.6:
                rows.append((time,np.arctanh(position_derivative/time_derivative)))
        write_table(data,table_name,{'t':[row[0] for row in rows],'eta':[row[1] for row in rows]})
    penrose_samples = sorted(set([math.tan(-math.pi/2+math.pi*index/2400) for index in range(1,2400)]
                                 +[-1e5,1e5,-1,0,1]))
    for name in CURVES:
        samples = penrose_samples
        if name == 'elliptic':
            samples = sorted(set(samples+[-peak_coordinate,peak_coordinate]))
        rows = []
        for coordinate in samples:
            time,position,time_derivative,position_derivative = worldline(name,coordinate,excursion)
            assert time_derivative > abs(position_derivative)
            x_value,t_value = compactify(time,position)
            assert abs(x_value)+abs(t_value) <= 1+1e-14
            assert not rows or t_value >= rows[-1][1]-1e-14
            rows.append((x_value,t_value))
        assert abs(rows[0][0]) < 1e-7 and abs(rows[0][1]+1) < 1e-4
        assert abs(rows[-1][0]) < 1e-7 and abs(rows[-1][1]-1) < 1e-4
        rows = [(0,-1)]+rows+[(0,1)]
        (data/('penrose_'+name+'.dat')).write_text('# X T\n'+''.join(
            f'{x_value:.15g} {t_value:.15g}\n' for x_value,t_value in rows),encoding='utf-8',newline='')
    print(f'Checked V={PEAK_SPEED}; elliptic L/tau={excursion:.12g}, t*/tau={peak_time:.12g}.')


def build_figures(work,engine):
    resolved = str(Path(engine).resolve()) if Path(engine).is_file() else shutil.which(engine)
    if not resolved:
        raise FileNotFoundError('pdfLaTeX was not found; supply --tex-engine /path/to/pdflatex.')
    for name in ('rapidity','penrose','shape_radiation'):
        source = FIGURES/(name+'.tex')
        shutil.copyfile(source,work/source.name)
        result = subprocess.run([resolved,'--disable-write18','-interaction=nonstopmode',
                                 '-halt-on-error','-file-line-error',source.name],cwd=work,
                                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                text=True,errors='replace',timeout=120)
        (work/(name+'_compile.txt')).write_text(result.stdout,encoding='utf-8')
        if result.returncode:
            raise RuntimeError('Figure build failed; inspect '+str(work/(name+'_compile.txt')))
        shutil.copyfile(work/(name+'.pdf'),FIGURES/(name+'.pdf'))
        print('Built figures/'+name+'.pdf')


def build_proof(work,engine,bibtex_engine=None):
    """Build the local proof in private scratch, keeping only its PDF publicly."""
    resolved = str(Path(engine).resolve()) if Path(engine).is_file() else shutil.which(engine)
    if not resolved:
        raise FileNotFoundError('pdfLaTeX was not found; supply --tex-engine /path/to/pdflatex.')
    source_folder = ROOT/'proof'
    proof_work = work/'proof'
    proof_work.mkdir(exist_ok=True)
    for path in source_folder.iterdir():
        if path.suffix in {'.tex','.bib'}:
            shutil.copyfile(path,proof_work/path.name)
    if not (proof_work/'proof.tex').is_file():
        raise FileNotFoundError('Missing proof/proof.tex.')
    transcripts = []

    def execute(arguments):
        result = subprocess.run(arguments,cwd=proof_work,stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,text=True,errors='replace',timeout=180)
        transcripts.append(result.stdout)
        (proof_work/'compile.txt').write_text('\n'.join(transcripts),encoding='utf-8')
        if result.returncode:
            raise RuntimeError('Proof build failed; inspect '+str(proof_work/'compile.txt'))

    command = [resolved,'--disable-write18','-interaction=nonstopmode',
               '-halt-on-error','-file-line-error','proof.tex']
    execute(command)
    if '\\bibdata{' in (proof_work/'proof.aux').read_text(encoding='utf-8',errors='replace'):
        if bibtex_engine is None:
            sibling = Path(resolved).with_name('bibtex'+Path(resolved).suffix)
            bibtex_engine = str(sibling) if sibling.is_file() else 'bibtex'
        bibtex = (str(Path(bibtex_engine).resolve()) if Path(bibtex_engine).is_file()
                  else shutil.which(bibtex_engine))
        if not bibtex:
            raise FileNotFoundError('BibTeX was not found; supply --bibtex-engine /path/to/bibtex.')
        execute([bibtex,'proof'])
    execute(command)
    execute(command)
    log = (proof_work/'proof.log').read_text(encoding='utf-8',errors='replace')
    if 'undefined references' in log or 'undefined citations' in log:
        raise RuntimeError('Unresolved proof references/citations; inspect '+str(proof_work/'proof.log'))
    shutil.copyfile(proof_work/'proof.pdf',source_folder/'proof.pdf')
    print('Built proof/proof.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tex-engine',default='pdflatex',help='pdfLaTeX executable name or full path.')
    parser.add_argument('--bibtex-engine',help='BibTeX executable name or full path.')
    parser.add_argument('--work-dir',type=Path,help='Keep generated tables and build files here.')
    parser.add_argument('--data-only',action='store_true',help='Generate/check samples without compiling.')
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--figures-only',action='store_true',help='Rebuild the three figure PDFs only.')
    selection.add_argument('--proof-only',action='store_true',help='Rebuild the theorem PDF only.')
    args = parser.parse_args()
    if args.proof_only and args.data_only:
        parser.error('--proof-only and --data-only cannot be combined.')
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='mirror-figures-'))
    work.mkdir(parents=True,exist_ok=True)
    if not args.proof_only:
        make_data(work)
    if not args.data_only:
        if not args.proof_only:
            build_figures(work,args.tex_engine)
        if not args.figures_only:
            build_proof(work,args.tex_engine,args.bibtex_engine)
    print('Working files: '+str(work))
