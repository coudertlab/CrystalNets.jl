"""Check all original-input equivalence certificates with Python's stdlib only.

No Julia, CrystalNets, Systre, NetworkX, SymPy, floats, or graph matching are used.
The certificate maps one connected periodic component type on each side.
"""
from collections import Counter,defaultdict
from itertools import combinations
from math import gcd
import hashlib,json,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent

def read_graph(path):
    values=list(map(int,path.read_text().split()));assert values[0]==3 and (len(values)-1)%5==0
    edges=[[values[i]-1,values[i+1]-1,*values[i+2:i+5]] for i in range(1,len(values),5)]
    return dict(nsites=max(max(e[:2]) for e in edges)+1,edges=edges)

def read_key(key):
    values=list(map(int,key.split()));d=values[0];assert d in (1,2)
    edges=[[values[i]-1,values[i+1]-1,*values[i+2:i+2+d],*([0]*(3-d))] for i in range(1,len(values),d+2)]
    return dict(nsites=max(max(e[:2]) for e in edges)+1,edges=edges)

def directed(g):
    edges=Counter()
    for u,v,*t in g['edges']:
        assert u!=v or any(t)
        edges[u,v,tuple(t)]+=1;edges[v,u,tuple(-x for x in t)]+=1
    assert all(v==1 for v in edges.values())
    return edges

def matvec(m,t):return tuple(sum(a*b for a,b in zip(row,t)) for row in m)
def det(m):
    if len(m)==1:return m[0][0]
    assert len(m)==2
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]

def lattice_divisor(columns,d):
    result=0
    if d==1:
        for v in columns:
            for x in v:result=gcd(result,x)
    else:
        for a,b in combinations(columns,2):
            for i,j in combinations(range(3),2):result=gcd(result,a[i]*b[j]-a[j]*b[i])
    return abs(result)

def coordinates(basis,v,d):
    if d==1:
        i=next(i for i in range(3) if basis[i][0]);numerator=v[i];denominator=basis[i][0]
        assert numerator%denominator==0;x=(numerator//denominator,)
    else:
        i,j=next((i,j) for i,j in combinations(range(3),2) if basis[i][0]*basis[j][1]-basis[i][1]*basis[j][0])
        denominator=basis[i][0]*basis[j][1]-basis[i][1]*basis[j][0]
        numerators=(v[i]*basis[j][1]-v[j]*basis[i][1],basis[i][0]*v[j]-basis[j][0]*v[i])
        assert all(a%denominator==0 for a in numerators)
        x=tuple(a//denominator for a in numerators)
    assert matvec(basis,x)==v
    return x

def intrinsic(g,certificate,d):
    q=certificate['tree_potential'];basis=certificate['basis'];n=g['nsites']
    assert len(q)==n and tuple(q[0])==(0,0,0) and len(basis)==3 and all(len(row)==d for row in basis)
    zero_adj=defaultdict(list);cycles=[];edges=[]
    for u,v,*t in g['edges']:
        cycle=tuple(t[i]+q[u][i]-q[v][i] for i in range(3));cycles.append(cycle)
        if not any(cycle):zero_adj[u].append(v);zero_adj[v].append(u)
        x=coordinates(basis,cycle,d)
        edges.append([u,v,*x,*([0]*(3-d))])
    reached={0};queue=[0]
    while queue:
        for v in zero_adj[queue.pop()]:
            if v not in reached:reached.add(v);queue.append(v)
    assert len(reached)==n,'The stated tree representatives must be reachable from the root'
    # Cycle vectors lie in the supplied lattice (checked above). Equality of
    # their determinantal divisors proves that this is the whole cycle lattice,
    # rather than a larger lattice that would introduce extra components.
    columns=[tuple(row[i] for row in basis) for i in range(d)]
    divisor=lattice_divisor(columns,d)
    assert divisor>0 and lattice_divisor(cycles,d)==divisor
    return dict(nsites=n,edges=edges)

def gauge(g,offsets,d):
    assert len(offsets)==g['nsites'] and all(len(v)==d for v in offsets)
    return dict(nsites=g['nsites'],edges=[[u,v,*(t[i]+offsets[v][i]-offsets[u][i] for i in range(d)),*([0]*(3-d))] for u,v,*t in g['edges']])

def cover(g,c,d):
    assert c['dimension']==d
    a,b,k=c['a'],c['b'],c['c'];assert a>0 and k>0 and 0<=b<a
    if d==1:assert b==0 and k==1
    n=g['nsites'];edges=[]
    for cell,(x,y) in enumerate((x,y) for y in range(k) for x in range(a)):
        for u,v,*t in g['edges']:
            dy=y+(t[1] if d==2 else 0);qy=dy//k;ry=dy%k
            dx=x+t[0]-b*qy;qx=dx//a;rx=dx%a
            edges.append([u+cell*n,v+(ry*a+rx)*n,qx,*([qy] if d==2 else []),*([0]*(3-d))])
    out=dict(nsites=n*a*k,edges=edges)
    assert len(directed(out))==a*k*len(directed(g))
    return out

def permutation_power(p,v,power):
    orbit=[v];w=p[v]
    while w!=v:orbit.append(w);w=p[w]
    return orbit[power%len(orbit)]

def verify(c,inputs):
    raw=[]
    for side in ('source','systre'):
        g=read_graph(inputs/c[side+'_file'])
        h=hashlib.sha256(json.dumps(g,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        assert h==c[side+'_hash'];raw.append(g)
    d=c['dimension'];assert d in (1,2)
    a=intrinsic(raw[0],c['source_intrinsic'],d);a=gauge(a,c['source_gauge'],d)
    b=intrinsic(raw[1],c['systre_intrinsic'],d);b=cover(b,c['systre_cover'],d);b=gauge(b,c['systre_cover_gauge'],d)
    n=a['nsites'];assert b['nsites']==n
    original_b=directed(b);generators=c['twist_generators'];assert len(generators)==d
    for p in generators:
        assert sorted(p)==list(range(n))
        assert Counter({(p[u],p[v],t):count for (u,v,t),count in original_b.items()})==original_b
    assert all(p[q[u]]==q[p[u]] for p in generators for q in generators for u in range(n))
    twisted=Counter()
    for (u,v,t),count in original_b.items():
        for i,p in enumerate(generators):v=permutation_power(p,v,-t[i])
        twisted[u,v,t]+=count
    assert all(twisted[v,u,tuple(-x for x in t)]==count for (u,v,t),count in twisted.items())
    witness=c['witness'];m=witness['basis'];mapping=witness['mapping'];r=witness['offsets']
    assert len(m)==d and all(len(row)==d for row in m) and abs(det(m))==1
    assert sorted(mapping)==list(range(n)) and len(r)==n and all(len(v)==d for v in r)
    mapped=Counter()
    for (u,v,t),count in directed(a).items():
        mt=matvec(m,t[:d]);ofs=tuple(mt[i]+r[v][i]-r[u][i] for i in range(d))+(0,)*(3-d)
        mapped[mapping[u],mapping[v],ofs]+=count
    assert mapped==twisted
    assert sum(mapped.values())==witness['edges_verified']
    return n,sum(mapped.values())

def verify_reduced(c):
    a,b=c['source'],c['target'];d=len(c['basis']);n=a['nsites']
    assert d in (1,2) and b['nsites']==n
    m=c['basis'];mapping=c['mapping'];r=c['offsets'];generators=c['twist_generators']
    assert len(m)==d and all(len(row)==d for row in m) and abs(det(m))==1
    assert sorted(mapping)==list(range(n)) and len(r)==n and all(len(v)==d for v in r)
    original_b=directed(b)
    assert len(generators)==d
    for perm in generators:
        assert sorted(perm)==list(range(n))
        assert Counter({(perm[u],perm[v],t):count for (u,v,t),count in original_b.items()})==original_b
    assert all(p[q[u]]==q[p[u]] for p in generators for q in generators for u in range(n))
    twisted=Counter()
    for (u,v,t),count in original_b.items():
        for i,perm in enumerate(generators):v=permutation_power(perm,v,-t[i])
        twisted[u,v,t]+=count
    assert all(twisted[v,u,tuple(-x for x in t)]==count for (u,v,t),count in twisted.items())
    mapped=Counter()
    for (u,v,t),count in directed(a).items():
        mt=matvec(m,t[:d]);ofs=tuple(mt[i]+r[v][i]-r[u][i] for i in range(d))+(0,)*(3-d)
        mapped[mapping[u],mapping[v],ofs]+=count
    assert mapped==twisted

def main():
    inputs=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent
    certificates=HERE/'input-certificates.json'
    rows=json.loads(certificates.read_text());assert [c['case_id'] for c in rows]==[f'D{i:02}' for i in range(1,24)]
    vertices=edges=0
    for c in rows:
        n,e=verify(c,inputs);vertices+=n;edges+=e
        print(c['case_id'],f'PASS: {n} vertices, {e} directed edges')
    print(f'PASS: all {len(rows)} original input pairs; {vertices} mapped vertices and {edges} directed edges.')
    cross=certificates.parent/'cross-case-certificate.json'
    if cross.exists():
        c=json.loads(cross.read_text());verify(c,inputs)
        assert c['source_systre_key']!=c['target_systre_key']
        byid={row['case_id']:row for row in rows}
        for case_id,key_field in [('D13','source_systre_key'),('D14','target_systre_key')]:
            assert directed(read_graph(inputs/byid[case_id]['systre_file']))==directed(read_key(c[key_field]))
        print('PASS: D13 and D14 source graphs are isomorphic despite distinct recorded Systre keys.')
    reductions=[line.split('\t') for line in (inputs/'reductions.tsv').read_text().splitlines() if not line.startswith('#')]
    reduced=json.loads((HERE/'reduced-certificates.json').read_text())
    assert [r['case_id'] for r in reduced]==[r[0] for r in reductions]
    for row,(_,example) in zip(reduced,reductions):
        certificate=row['certificate'];verify_reduced(certificate)
        pair=[directed(read_key((inputs/(example+'_'+side+'.pgraph')).read_text())) for side in ('a','b')]
        certified=[directed(certificate[side]) for side in ('source','target')]
        assert (pair==certified or pair==certified[::-1]),'Certificate must match the shipped reduced graphs'
    print('PASS: all 23 reduced-pair certificates match the shipped fixtures.')
if __name__=='__main__':main()
