"""Semantic adversaries: consistent arithmetic with invalid mathematical contracts."""
import copy,gzip,json,os
from pathlib import Path
import verify_flint as v

def reject(label,fn,accepted):
    try:fn()
    except (AssertionError,ValueError,ZeroDivisionError):accepted.append(label)
    else:raise AssertionError('semantic mutant accepted: '+label)

def gram(rec,H,sm):
    n=rec['degree'];xs=list(map(v.Q,rec['x_grid']));r=rec['channels'][sm];K=r['count'];V=r['vectors'];A=[[0]*K for _ in range(K)];m=int(sm);k=abs(n-m);j=n+m;x=xs[0];g=v.G(n,H,x)
    org=[v.wholeup(k+3*x*g/(2*(k+1)),v.D),v.wholeup(j+3*x*g/(2*(j+1)),v.D)]
    for p in range(K):
        for q in range(K):A[p][q]+=org[0]*V[0][p]*V[0][q]+org[1]*V[1][p]*V[1][q]
    for i,(a,b) in enumerate(zip(xs[:-1],xs[1:])):
        d,c=v.lowercell(n,H,m,a,b);M=[[d[0],d[1],c[0],c[1]],[d[1],d[2],c[1],c[2]],[c[0],c[1],d[0],d[1]],[c[1],c[2],d[1],d[2]]];w=V[2*i:2*i+4]
        for p in range(K):
            for q in range(K):A[p][q]+=sum(M[s][t]*w[s][p]*w[t][q] for s in range(4) for t in range(4))
    r['gram_integer']=[[str(x) for x in row] for row in A];r['leading_minors']=[str(v.fmpz_mat([row[:sz] for row in A[:sz]]).det()) for sz in range(1,K+1)]

def main():
    root=Path(__file__).parent;accepted=[]
    for n in v.EXPECTED_INDICES:
        profile=json.loads((root/f'profiles/profile_n{n}.json').read_text());hs=v.profiles(profile)
        lo=json.loads((root/f'certificates/lower_n{n}.json').read_text());up=json.loads(gzip.decompress(next((root/'certificates').glob(f'upper_n{n}*.gz')).read_bytes()))
        bad=copy.deepcopy(lo)
        for sm in bad['channels']:gram(bad,hs['lower'],sm)
        reject(f'n{n}: wrong-side lower with consistently recomputed Gram',lambda:v.lower(bad,hs['lower']),accepted)
        reject(f'n{n}: wrong-side upper profile',lambda:v.upper(up,hs['upper']),accepted)
        lc={m:r['count'] for m,r in lo['channels'].items()};uc={m:r['count'] for m,r in up['channels'].items()}
        reject(f'n{n}: unmatched paired counts',lambda:v.matching_counts(n,lo['lower_index'],up['upper_index']+2,lc,uc),accepted)
        wrong=copy.deepcopy(up);wrong['channels']['2']['tail']['robin']=40
        reject(f'n{n}: unjustified exterior Robin constant',lambda:v.upper(wrong,hs['lower']),accepted)
        wrong=copy.deepcopy(lo);r=wrong['channels']['2'];r['count']+=1
        for row in r['vectors']:row.append(row[0])
        gram(wrong,hs['upper'],'2');wrong['lower_index']+=2
        reject(f'n{n}: dependent extra trial vector with consistent Gram',lambda:v.lower(wrong,hs['upper']),accepted)
        wrong=copy.deepcopy(up);wrong['x_grid'][-1]=str(v.Q(wrong['x_grid'][-1])-1)
        reject(f'n{n}: stale last-cell geometry',lambda:v.upper(wrong,hs['lower']),accepted)
    n=5;profile=json.loads((root/'profiles/profile_n5.json').read_text());hs=v.profiles(profile);y=v.P([0,1])
    for side,delta in [('lower',v.Q(1,1000)),('upper',v.Q(-1,1000))]:
        H=hs[side]+delta*(1-y)**2
        reject(f'n5: {side} shifted toward profile with residual recomputed',lambda H=H,side=side:v.checked_barrier(5,H,side),accepted)
    out=Path(os.environ.get('SM_OUT','.'));(out/'semantic_mutants.json').write_text(json.dumps({'rejected':accepted,'count':len(accepted)},indent=2));print('SEMANTIC MUTANTS PASS',len(accepted),flush=True)
if __name__=='__main__':main()
