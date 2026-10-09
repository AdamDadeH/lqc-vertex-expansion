"""Minimal parser for Mathematica .nb files -> readable Markdown digest.

Renders Input cells as (roughly InputForm) Mathematica code, Text/Title/Section cells as
Markdown, and optionally Output/Print cells as truncated comments.  No Mathematica needed.

Standalone: python archive/nbdigest.py archive/mathematica --out -o archive/digests
"""
import re, sys, os, collections

class Sym(str): pass

class Call:
    def __init__(self, head, args): self.head, self.args = head, args
    def __repr__(self): return f"{self.head}[{', '.join(map(repr,self.args))}]"

class P:
    def __init__(self, s): self.s=s; self.i=0; self.n=len(s)
    def ws(self):
        while True:
            while self.i<self.n and self.s[self.i].isspace(): self.i+=1
            if self.s.startswith("(*", self.i):
                j=self.s.find("*)", self.i+2); self.i = j+2 if j>=0 else self.n
            elif self.s.startswith("\\\n", self.i) or self.s.startswith("\\\r\n", self.i):
                self.i+=2
            else: return
    def peek(self):
        self.ws(); return self.s[self.i] if self.i<self.n else ''
    def expect(self,c):
        self.ws()
        if self.s[self.i]!=c: raise SyntaxError(f"expected {c!r} at {self.i}: {self.s[self.i:self.i+60]!r}")
        self.i+=1
    def string(self):
        assert self.s[self.i]=='"'; self.i+=1
        out=[]
        while True:
            c=self.s[self.i]
            if c=='"': self.i+=1; return ''.join(out)
            if c=='\\':
                d=self.s[self.i+1]
                if d=='[':
                    j=self.s.index(']',self.i); out.append(self.s[self.i:j+1]); self.i=j+1; continue
                if d=='n': out.append('\n'); self.i+=2; continue
                if d=='t': out.append('\t'); self.i+=2; continue
                if d=='\n': self.i+=2; continue
                if d=='<' or d=='>': out.append(''); self.i+=2; continue
                if d=='!' or d=='(' or d==')' or d=='@' or d=='*' or d=='_' or d=='^' or d=='%' or d=='&' or d=='+': out.append('\\'+d); self.i+=2; continue
                out.append(d); self.i+=2; continue
            out.append(c); self.i+=1
    def atom(self):
        self.ws()
        c=self.s[self.i]
        if c=='"': return self.string()
        if c=='{':
            self.i+=1; items=[]
            if self.peek()=='}': self.i+=1; return items
            while True:
                items.append(self.expr())
                if self.peek()==',': self.i+=1; continue
                self.expect('}'); return items
        if c=='(':
            self.i+=1; e=self.expr(); self.expect(')'); return e
        m=re.compile(r'[-+]?(\d+\.?\d*|\.\d+)(`[\d.]*)?(\*\^[-+]?\d+)?|\d+\^\^[0-9a-zA-Z.]+').match(self.s,self.i)
        if m and not (c in '-+' and not self.s[self.i+1:self.i+2].strip()):
            self.i=m.end(); return m.group(0)
        m=re.compile(r'([A-Za-z$]|\\\[\w+\])([A-Za-z0-9$]|\\\[\w+\])*(`([A-Za-z0-9$`]|\\\[\w+\])*)?').match(self.s,self.i)
        if m:
            self.i=m.end(); name=Sym(m.group(0))
            # call(s)
            e=name
            while self.i<self.n and self.s[self.i]=='[':
                self.i+=1; args=[]
                if self.peek()==']': self.i+=1; e=Call(e,args); continue
                while True:
                    args.append(self.expr())
                    if self.peek()==',': self.i+=1; continue
                    self.expect(']'); break
                e=Call(e,args)
            return e
        if c=='-':
            self.i+=1; return Call(Sym('Minus'),[self.atom()])
        if c=='#': self.i+=1; return Sym('#')
        raise SyntaxError(f"bad atom at {self.i}: {self.s[self.i:self.i+80]!r}")
    def expr(self):
        e=self.atom()
        while True:
            p=self.peek()
            if self.s.startswith('->',self.i) or self.s.startswith(':>',self.i):
                op=self.s[self.i:self.i+2]; self.i+=2; e=Call(Sym(op),[e,self.expr()])
            elif self.s.startswith('&',self.i):
                self.i+=1; e=Call(Sym('Function'),[e])
            elif p in '*/+-^' and p:
                self.i+=1; e=Call(Sym(p),[e,self.atom()])
            elif p and (p.isalnum() or p in '$({"\\'):
                e=Call(Sym('Times'),[e,self.atom()])   # implicit multiplication
            else: return e

NAMED = {'Theta':'θ','Delta':'δ','Pi':'π','ImaginaryI':'I','IndentingNewLine':'\n','Rule':'->','RuleDelayed':':>',
 'Infinity':'Infinity','Phi':'φ','Psi':'ψ','Lambda':'λ','Mu':'μ','Nu':'ν','Alpha':'α','Beta':'β','Gamma':'γ','Epsilon':'ε',
 'Sigma':'σ','Omega':'ω','Rho':'ρ','Tau':'τ','Xi':'ξ','Eta':'η','Kappa':'κ','Zeta':'ζ','Chi':'χ','CapitalDelta':'Δ','CapitalPhi':'Φ',
 'CapitalPsi':'Ψ','CapitalTheta':'Θ','CapitalGamma':'Γ','CapitalLambda':'Λ','CapitalSigma':'Σ','CapitalOmega':'Ω','CapitalPi':'Π',
 'Equal':'==','NotEqual':'!=','LessEqual':'<=','GreaterEqual':'>=','Element':'∈','ExponentialE':'E','Degree':'Degree',
 'InvisibleSpace':'','InvisibleApplication':'','Times':'*','LeftDoubleBracket':'[[','RightDoubleBracket':']]','Conjugate':'^*',
 'CurlyPhi':'ϕ','CurlyEpsilon':'ϵ','CurlyTheta':'ϑ','DifferentialD':'d','PartialD':'D','Sum':'Σ','Product':'Π','Integral':'∫',
 'Dash':'-','LongEqual':'==','RawWedge':'^','Wolfram':'','Mho':'℧','Bullet':'•','LeftArrow':'<-','RightArrow':'->','Cross':'×',
 'Transpose':'^T','HBar':'ħ','Null':'','Prime':"'",'DoubleStruckCapitalR':'ℝ','Ellipsis':'...','ScriptCapitalH':'ℋ','ScriptL':'ℓ',
 'SpanFromLeft':'','Dagger':'†','Backslash':'\\','LineSeparator':'\n','NewLine':'\n','Proportional':'∝','Tilde':'~','Not':'!','And':'&&','Or':'||','Implies':'=>','TildeTilde':'≈','UndirectedEdge':'<->','DirectedEdge':'->'}
def named(s):
    return re.sub(r'\\\[(\w+)\]', lambda m: NAMED.get(m.group(1), '\\['+m.group(1)+']'), s)

def box(b):
    """Render a box expression to text (roughly InputForm)."""
    if isinstance(b,str) and not isinstance(b,Sym): return named(b)
    if isinstance(b,Sym): return str(b)
    if isinstance(b,list): return ''.join(box(x) for x in b)
    h=str(b.head); a=b.args
    if h=='RowBox': return box(a[0])
    if h=='FractionBox': return f"({box(a[0])})/({box(a[1])})"
    if h=='SuperscriptBox': return f"({box(a[0])})^({box(a[1])})"
    if h=='SubscriptBox': return f"Subscript[{box(a[0])}, {box(a[1])}]"
    if h=='SubsuperscriptBox': return f"Subscript[{box(a[0])}, {box(a[1])}]^({box(a[2])})"
    if h=='SqrtBox': return f"Sqrt[{box(a[0])}]"
    if h=='RadicalBox': return f"({box(a[0])})^(1/({box(a[1])}))"
    if h=='UnderoverscriptBox': return f"{box(a[0])}_({box(a[1])})^({box(a[2])})"
    if h=='UnderscriptBox': return f"{box(a[0])}_({box(a[1])})"
    if h=='OverscriptBox': return f"Overscript[{box(a[0])}, {box(a[1])}]"
    if h=='GridBox':
        rows=a[0]; return "{" + ", ".join("{"+", ".join(box(c) for c in r)+"}" for r in rows) + "}"
    if h in ('StyleBox','TagBox','InterpretationBox','FormBox','AdjustmentBox','ItemBox','PaneBox','FrameBox','TooltipBox','ButtonBox','DynamicBox','TemplateBox','PanelBox'): return box(a[0])
    if h=='GraphicsBox': return "<<Graphics>>"
    if h=='Graphics3DBox': return "<<Graphics3D>>"
    if h=='ErrorBox': return "<<Error>>"
    if h=='TextData': return box(a[0])
    if h=='Cell': return box(a[0])
    if h=='BoxData': return box(a[0])
    if h=='Minus': return '-'+box(a[0])
    return f"<<{h}>>"

def cells(node, out):
    if isinstance(node,Call):
        if str(node.head)=='Cell' and node.args:
            content=node.args[0]; style=node.args[1] if len(node.args)>1 and isinstance(node.args[1],str) and not isinstance(node.args[1],Sym) else '?'
            if isinstance(content,Call) and str(content.head)=='CellGroupData':
                for c in content.args[0]: cells(c,out)
                return
            out.append((style, content)); return
        for x in node.args: cells(x,out)
    elif isinstance(node,list):
        for x in node: cells(x,out)

def digest(path, include_output=False, max_out=400):
    s=open(path,encoding='utf-8',errors='replace').read()
    s=s[s.index('Notebook['):]
    nb=P(s).expr()
    out=[]; cells(nb,out)
    lines=[]
    counts=collections.Counter(st for st,_ in out)
    for st,content in out:
        txt=box(content)
        if st in ('Title','Section','Subsection','Subsubsection','Chapter','Subtitle'):
            lines.append(f"\n{'#'*{'Title':1,'Chapter':1,'Subtitle':2,'Section':2,'Subsection':3,'Subsubsection':4}[st]} {txt.strip()}\n")
        elif st=='Text': lines.append(f"> {txt.strip()}\n")
        elif st=='Input': lines.append("```mathematica\n"+txt.strip()+"\n```")
        elif st=='Output' and include_output:
            t=txt.strip(); 
            if len(t)>max_out: t=t[:max_out]+f" ...[{len(txt)} chars]"
            lines.append("(* Out: "+t+" *)")
        elif st in ('Print','Message') and include_output: lines.append("(* "+st+": "+txt.strip()[:200]+" *)")
    return counts, '\n'.join(lines)


def main(argv=None):
    import argparse, pathlib
    ap=argparse.ArgumentParser(description="Render Mathematica notebooks as Markdown digests.")
    ap.add_argument('paths', nargs='+', help='.nb files or directories (searched recursively)')
    ap.add_argument('--out', action='store_true', help='include Output/Print cells (truncated)')
    ap.add_argument('--max-out', type=int, default=300, help='truncate outputs to this many chars')
    ap.add_argument('-o', '--outdir', help='write one .md per notebook here instead of stdout')
    a=ap.parse_args(argv)
    files=[]
    for p in a.paths:
        p=pathlib.Path(p)
        files += sorted(p.rglob('*.nb')) if p.is_dir() else [p]
    for f in files:
        counts,text=digest(str(f), a.out, a.max_out)
        if a.outdir:
            out=pathlib.Path(a.outdir)/(f.name[:-3]+'.md')
            out.write_text(f"# {f.name}\n\n_cells: {dict(counts)}_\n\n"+text+"\n")
            print(f"wrote {out}  {dict(counts)}")
        else:
            print(f"==== {f}  {dict(counts)}"); print(text)

if __name__=='__main__':
    main()
