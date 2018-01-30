with(combinat):with(ListTools): with(ArrayTools): with(LinearAlgebra):with(convex): 

# For a given monomial, the number of mons with degree difference being at most 1 is:
# sum_{k=1}^floor(n/2) (2*k)! choose(n,2k) + 2 * sum_{k=1}^ceil(n/2) (2*k-1)! choose(n,2k-1) 
# something close to n/2 * factorial(n)

monshom := proc(n,d)
  local l1;
  l1 := [seq(1,i=1..n)];
  return map(_c -> _c - l1, composition(n+d, n)):
end;

# For both revmons functions, the bottlneck is the composition function from package combinat

revmons := proc(n,d) 
  if d = 0 then 
  return [[seq(0,i=1..n)]]:
  else return [op(monshom(n,d)), op(revmons(n,d-1))]:
  fi:
end;

revmons2 := proc(n,d) 
  local mons0;
  mons0 := [[seq(0,i=1..n)]]:
  return [seq(op(monshom(n,d-i+1)),i=1..d),mons0]
end;

monspt := proc(n,d)
  return Reverse(revmons(n,d)):
end;

mons := proc(n,d,X,ms::list := monspt(n,d))
  return map(a->mul(X[i]^a[i],i=1..n),ms);
#  return [seq(mul(X^~a), a in ms)];
end;

benchmons := proc(ms)
  local m, i;
  #ms := mons(n,d);
  for i from 1 to nops(ms) do
    m := ms[i];
    Search(m,ms):
  od: 
  return 0:
end;

support := proc (f,X,monslist)
  local cf, idxf,monsf;
  cf := coeffs(f,X,'monsf');
  monsf := [monsf];
  idxf := seq(Search(m,monslist),m in monsf);
  return [cf],[idxf];
end;

dense_perturbation := proc(ms)
  return add(m^2,m in ms); 
end;

sos2sdp := proc(r,X,relaxorder, mspt,msk,rc,ridx,precSVD::integer:=10,precSDP::integer:=200,epsStar::integer:=3,precRound::integer:=30,dataprec::integer:=10,gmp::boolean:=false,algo::integer:=1,g::list:=[],gc::list:=[],gidx::list:=[])
  local start_idx, sub_idx,nvars,zero,d,k,n,ng,r0,msdp,nsdp,nblock,fd,rList,i,rfloat,j,lowerbnd,eigs,eivs,eigs2,eivs2,Y, Yrat, absorb_tbl,m, absorb_PP,nsdpg,nsdponeg,gic,giidx,ni,ig,igc,gifloat,mgi,eigsi,eivsi;
  k := relaxorder; d = 2*k: n := nops(X); ng := nops(g): zero := [seq(0,i=1..n)]:
  msdp := nops(mspt); nsdp := nops(msk); nblock := 1+ng: 
  fd := fopen("in.dat-s",WRITE,TEXT); 
  rList := Array([seq(0,i=1..msdp)]):
  for i from 1 to nops(rc) do rList[ridx[i]]:=rc[i]: od:

  if mspt[1] = zero then nvars := msdp - 1: start_idx := 2: r0 := rList[1]: 
 else sub_idx := Search(2*mspt[1],mspt)-1: nvars := msdp - sub_idx: start_idx := 1: r0 := 0: fi: # test whether 0 belongs the multi-index list
  writeline(fd,convert(nvars,string)); # Number of moment variables - 1 (since y0 = 1)
  writeline(fd,convert(nblock,string)); # Number of blocks: 1 for the moment matrix
  fprintf(fd,"%d ",nsdp): # Size of the first block: the moment matrix has size nsdp
  nsdpg := locmatsizes(g,k,n):
  #lprint(nsdpg);
  for ni in nsdpg do fprintf(fd,"%d ",ni): od: # Size of the localizing matrices: for a polynomial gi of degree di, the size is binomial(n+k-floor(di/2),n)
  fprintf(fd,"\n");

  rfloat := convert(rList,float,dataprec):
  #lprint(rfloat);

  if mspt[1] = zero then 
    for i from start_idx to msdp do fprintf(fd,"%s ",convert(rfloat[i],string));od;
  else
    for i from sub_idx + 1 to msdp do fprintf(fd,"%s ",convert(rfloat[i],string));od;
  fi:
  fprintf(fd,"\n");

### The case when y0 = 1 because 0 belongs to the multi-index list
# starts with the first row of the moment matrix
  if mspt[1] = zero then 
  writeline(fd,"0 1 1 1 -1"): # Moment variable y0 = 1
  for j from 2 to nsdp do
    fprintf(fd, "%d 1 %d %d 1\n", j-1, 1, j); # Moment variable y_{j-1} at first row and column j
  od;
# then the other rows of the moment matrix
  for i from 2 to nsdp do
    for j from i to nsdp do
      k := Search(mspt[i]+mspt[j],mspt):
      fprintf(fd, "%d 1 %d %d 1\n", k-1, i, j); # Moment variable y_{k-1} at row i and column j
    od:
  od:
# then the localizing matrices
for ig from 1 to ng do
   gic := gc[ig]: giidx := gidx[ig]: 
   for igc from 1 to nops(gic) do
     mgi := mspt[giidx[igc]]:
     #lprint(mgi);
     gifloat := convert(gic[igc],float,dataprec):
     if mgi = zero then 
      fprintf(fd, "0 %d 1 1 %s\n", ig+1,convert(-gifloat,string)); # Moment variable y_{j-1} at first row and column j
      for j from 2 to nsdpg[ig] do
        fprintf(fd, "%d %d 1 %d %s\n", j-1, ig+1, j, convert(gifloat,string)); # Moment variable y_{j-1} at first row and column j
      od;
     for i from 2 to nsdpg[ig] do
       for j from i to nsdpg[ig] do
         k := Search(mspt[i]+mspt[j]+mgi,mspt):
         fprintf(fd, "%d %d %d %d %s\n",k-1, ig+1, i, j, convert(gifloat,string)); # Moment variable y_{k-1} at row i and column j
       od:
     od: 
     else
     for i from 1 to nsdpg[ig] do
       for j from i to nsdpg[ig] do
         k := Search(mspt[i]+mspt[j]+mgi,mspt):
         fprintf(fd, "%d %d %d %d %s\n",k-1, ig+1, i, j, convert(gifloat,string)); # Moment variable y_{k-1} at row i and column j
       od:
     od: 
     fi:
   od:
od:

### The case when y0 = 1 because 0 does not belong to the multi-index list
# there is a bug in sdpa-gmp when not subtracting the "useless" moment variables with k - sub_idx
  else
  for i from 1 to nsdp do
    for j from i to nsdp do
      k := Search(mspt[i]+mspt[j],mspt):
      fprintf(fd, "%d 1 %d %d 1\n", k-sub_idx, i, j); # Moment variable y_{k-sub_idx} at row i and column j
    od;
  od;
  fi:

#### Loop again to find the absorbers (this loop can be saved thanks to the previous one)
  absorb_tbl := table([seq(0,i=1..msdp)]):
  for i from 1 to nsdp do
    for j from i to nsdp do
      m := mspt[i]+mspt[j]:
            if not_even(m) then k:= Search(m,mspt): absorb_tbl[k]:=i: fi:
    od:
  od:

#### In this table, for each gamma (appearing at index k in mspt), one indicates all monomials alpha (and implicitely beta) such that gamma = alpha + beta 
  absorb_PP := table([seq([],i=1..msdp)]):
  for i from 1 to nsdp do
  k:= Search(2*mspt[i],mspt): absorb_PP[k]:=[op(absorb_PP[k]),i,i]:
    for j from i+1 to nsdp do
      m := mspt[i]+mspt[j]: k:= Search(m,mspt): absorb_PP[k]:=[op(absorb_PP[k]),i,j,j,i]:
    od:
  od:
  #lprint(msdp); lprint(mspt); lprint(absorb_PP);

#### old version, works only when all monomials belong to the NP
#  for j from 2 to nsdp do
#    fprintf(fd, "%d 1 %d %d 1\n", j-1, 1, j); # Moment variable y_{j-1} at first row and column j
#  od;
#  for i from 2 to nsdp do
#    for j from i to nsdp do
#      k := Search(mspt[i]+mspt[j],mspt):
#      fprintf(fd, "%d 1 %d %d 1\n", k-1, i, j); # Moment variable y_{k-1} at row i and column j
#    od;
#  od;

  fclose(fd);
  system("sed -i 's/ \\./ 0\\./g' in.dat-s"); system("sed -i 's/^\\./0\\./g' in.dat-s"); system("sed -i 's/^-\\./-0\\./g' in.dat-s"); system("sed -i 's/ -\\./ -0\\./g' in.dat-s");
  system("rm -f out.dat-s"); system("rm -f out.mm");
  if not gmp then
    system("sdpa -ds in.dat-s -o out.dat-s -p param.sdpa > /dev/null");
  else
    write_param(precSDP,epsStar,epsStar);
    system("sdpa_gmp -ds in.dat-s -o out.dat-s -p my_param_gmp.sdpa > /dev/null");
  fi:
  system("echo $(grep objValPrimal out.dat-s) ';' 'yMat:=' $(sed -n '/yMat/,/main/{//!p}' out.dat-s) ';' >> out.mm");
  #system("echo 'yMat:=' $(sed -n '/yMat/,/main/{//!p}' out.dat-s) ';' >> out.mm");
  system("sed -i 's/ =/ :=/g' out.mm"): 
  system("sed -i 's/{/[/g' out.mm"): system("sed -i 's/}/]/g' out.mm"); system("sed -i 's/] \\[/], \\[/g' out.mm");
  read "out.mm": 
  if algo = 1 then lowerbnd := r0 + convert(objValPrimal,rational,exact) : else lowerbnd := 0: fi:
  printf("lower bound: "); lprint(evalf(lowerbnd));
  eigs:=Array([]); eivs := Array([]);
  nsdponeg:=[nsdp,op(nsdpg)]:
  for i from 1 to ng+1 do
    Y := Matrix(yMat[i]): 
    if algo = 1 then Yrat:=Matrix(Y): 
    else # case when algo = 2
     if precRound >= 30 then Y := convert(Y,rational,exact): else Y := convert(Y,rational,precRound) fi:
     Y := (Y + Transpose(Y))/2;
     if ng = 0 then Yrat := absorber_PP(Y,rList,absorb_PP,nsdp,mspt): else Yrat:=Y: fi:
     printf("err="); lprint(expand(Transpose(Vector(msk)).Yrat.Vector(msk)+lowerbnd-r)) 
    fi:
    (eigsi,eivsi) := eigseivs(Yrat,X,msk[1..nsdponeg[i]],precSVD,precRound,algo);
    eigs:=Concatenate(1,eigs,eigsi):  eivs:=Concatenate(1,eivs,eivsi):
   od:
  return convert(eigs,list), convert(eivs,list),lowerbnd,absorb_tbl,nsdponeg:
end;

absorber_PP := proc(Yrat,rList,absorb_PP,nsdp,mspt)
  local Yproj,i,j,i1,k,m,list_m,nopsm,maxY,minY,yij;
  # for i from 1 to nsdp do
  #  for j from i to nsdp do 
  #    yij := Yrat[i,j]: if abs(yij) < 1/10^3 then Yrat[i,j] := 0: Yrat[j,i]:=0: fi:
  #  od:
  #od:
   #Yproj := Matrix(Yrat,shape='symmetric');
   Yproj := Matrix(nsdp,nsdp,shape='symmetric');
   maxY := max(max(abs(Yproj))):    minY := min(min(abs(Yproj))): 
   #lprint(evalf(minY)); lprint(evalf(maxY)); #lprint(Y); #lprint(rList);
   for i from 1 to nsdp do
    for j from i to nsdp do 
      m := mspt[i] + mspt[j]; k:= Search(m,mspt):
      list_m := [op(absorb_PP[k])]: nopsm := nops(list_m)/2:
      Yproj[i,j] := Yrat[i,j] - 1/nopsm*(add(Yrat[list_m[2*im-1],list_m[2*im]],im=1..nopsm) - rList[k]):
    od:
  od:
  return Yproj: 
end:


write_param := proc(precSDP,epsStar,epsDash)
  local fd;
  fd := fopen("my_param_gmp.sdpa",WRITE,TEXT); 
  fprintf(fd,"300	unsigned int maxIteration;\n");
  fprintf(fd,"1.0E-%d	\t double 0.0 < epsilonStar;\n",epsStar);
  fprintf(fd,"1.0E5   double 0.0 < lambdaStar;\n");
  fprintf(fd,"2.0   	double 1.0 < omegaStar;\n");
  fprintf(fd,"-1.0E5  double lowerBound;\n");
  fprintf(fd,"1.0E5   double upperBound;\n");
  fprintf(fd,"0.1     double 0.0 <= betaStar <  1.0;\n");
  fprintf(fd,"0.3     double 0.0 <= betaBar  <  1.0, betaStar <= betaBar;\n");
  fprintf(fd,"0.9     double 0.0 < gammaStar  <  1.0;\n");
  fprintf(fd,"1.0E-%d	\t double 0.0 < epsilonDash;\n",epsDash);
  fprintf(fd,"%d     precision\n",precSDP);
  fclose(fd);
end;

checkrational := proc(U)
  local v:
  for v in U do:
    if not type(convert(v,rational),realcons) then 
    lprint(v): 
    error "Non Rational Cholesky factor, retry with gmp = true":fi:
  od:
  return:
end;

eigseivs := proc(Yrat,X,ms,precSVD,precRound,algo)
  local v, e, msvec, eigs, eivs, U,S,V, Ysvd,SVD,ti,tcmp;
  SVD := false:
  msvec := Vector(ms);
  # if SVD then 
  #  Digits:=precSVD;
  #  lprint("starting SVD");
  #  ti := time():
  #  (U,S,V) := MTM[svd](Yrat); 
  #  tcmp := time() - ti: lprint (tcmp);
  #  lprint("ending SVD");
  #else
    ti := time():
    lprint("starting Cholesky");
    if algo = 1 then Digits := precSVD:fi:
    U := LUDecomposition(Yrat,method='Cholesky');
    checkrational(U):
    tcmp := time() - ti:
    #lprint(U);
    lprint (tcmp);
    lprint("ending Cholesky");
    S := IdentityMatrix(nops(ms));
  # fi:
    eigs := Diagonal(S); eivs := Transpose(Transpose(msvec).U);
    Digits:=10;
    Ysvd := U.S.V^%T;
    #lprint(Ysvd); #lprint(max(max(Y),-min(Y))); #lprint(max(max(Ysvd),-min(Ysvd))); #lprint(max(Y - Ysvd));
    if precRound >= 30 then lprint("exact");  eigs := convert(eigs,rational,exact); eivs := map (_e -> convert(_e,rational,exact), eivs): else eigs := convert(eigs,rational,precRound); eivs := map (_e -> convert(_e,rational,precRound), eivs); fi:
    return (eigs, eivs):
end;


neighbours1 := proc(n)
  local i,nones,nmones,nzeros,nlist,n1;
  nlist := [];
  for i from 1 to floor(n/2) do
    nones:=i; nmones:=i; nzeros:=n-2*i; nlist:=[op(nlist),op(permute([1$nones, 0$nzeros, -1$nmones]))];
  od:
  for i from 1 to ceil(n/2) do
    nones:=i; nmones:=i-1; nzeros:=n-2*i+1; n1 := permute([1$nones, 0$nzeros, -1$nmones]);
    nlist := [op(nlist),op(n1),op(-n1)]:
  od:
  return nlist:
end;

count_n1 := proc(n)
  add(multinomial(n,i,i,n-2*i),i=1..floor(n/2)) + 2*add(multinomial(n,i,i-1,n-2*i+1),i=1..ceil(n/2));
end;

even_n1 := proc(n,mspt)
  local n1,a,list_n1,cardn1,n1a,b,ab;
  list_n1 := [];
  n1 := neighbours1(n); 
  cardn1 := nops(n1);
  for a in mspt do:
    n1a := [];
    for b from 1 to cardn1 do 
      ab := a+n1[b]: if not has(ab,-1) then n1a := [op(n1a),ab]:fi:
    od:
    list_n1:=[op(list_n1),n1a]:
  od:
  return list_n1:
end;

# The following degree_diff1 procedure returns wrong results: for instance Y^2 does not belong to neighbours1(X^2) but degree_diff1(X^2,Y^2) = true
# degree_diff1 := proc(a,b)
#  local c := b - a;
#  if abs(add(i, i in c))<= 1 then return true else return false:fi:
#end;

not_even := proc(m)
  return has(1,map(_c -> irem(_c,2), m)):
end;

### This procedure sometimes fails to find absorbing polynomials since the resulting decomposition yields monomials whose squares do not necessarily belong to the NP
decomp_mon := proc(a,n)
  local c,b,degb,cm,cnt,i;
  c := map(_c -> iquo(_c,2), a);
  b := a - 2*c;
  degb := ceil(add(i,i in b)/2);
  cm := [seq(0,i=1..n)]:
  cnt := 0;
  for i from 1 to n do:
    if b[i] = 1 then cm[i] := b[i]: cnt := cnt + 1: fi:
    if cnt >= degb then break: fi:
  od:
  return c,cm,b - cm;
end;

truncate_withNP := proc(f,X,mspt,msptk)
  local NP,msptkNP,msptNP,m,i;
  NP := newtonpolytope(f, X):
  msptkNP := [];
  for m in msptk do if contains(NP,2*m) then msptkNP := [op(msptkNP),m] : fi: od:
  msptNP := [op(msptkNP)]:
  for i from nops(msptk) + 1 to nops(mspt) do m := mspt[i]:  if contains(NP,m) then msptNP := [op(msptNP),m] : fi: od:
  return msptNP,msptkNP:
end:

testNP := proc(f,X)
  local d,k,n,card_nk,mspt,msptk;
  d := degree(f): k := d/2: n := nops(X): card_nk := binomial(n+k,k):
  mspt := monspt(n,d): msptk := mspt[1..card_nk]:
  return truncate_withNP(f,X,mspt,msptk):
end:

absorber := proc(u,X,e,even_mons,ms,mspt,absorb_tbl)
  local i,j,k,ucoeffs,uidx,uc,err_list,err,m,bad_m,n,c,cm,cp,cfs,sqs,n1,m1,m2,k1,k2;
  #printf("\nu = "); lprint(evalf(u));
  ucoeffs,uidx := support(u,X,ms);
  n := nops(X): 
  #printf("\neven_mons = "); lprint(even_mons);
  err_list := Array([seq(e,i=1..nops(even_mons))]); 
  cfs := []; sqs := [];
  for i from 1 to nops(ucoeffs) do
    uc := ucoeffs[i]; bad_m := mspt[uidx[i]]:
    if not_even(bad_m) then
      k := Search(bad_m,mspt):m1:=mspt[absorb_tbl[k]]:m2:=bad_m - m1;
      #printf("\n"); lprint(bad_m); lprint(m1); lprint(m2);
      cfs := [op(cfs),abs(uc)/2]; 
      sqs := [op(sqs),mul(X[i]^m1[i],i=1..n)+sign(uc)*mul(X[i]^m2[i],i=1..n)]:
      k1:=Search(2*m1,even_mons): k2:=Search(2*m2,even_mons):
      err_list[k1] := err_list[k1]-1/2*abs(uc): err_list[k2] := err_list[k2]-1/2*abs(uc):

    else 
      k := Search(bad_m,even_mons): err_list[k]:=err_list[k]+uc:
    fi:
  od:
  err_list := convert(err_list,list);
  # printf("\nerr list = "); lprint(evalf(err_list));
  cfs := [op(cfs),op(err_list)]; sqs := [op(sqs),seq(mul(X[i]^(m[i]/2),i=1..n), m in even_mons)];
  return cfs,sqs:
end;


old_absorber := proc(u,X,e,even_mons,ms,mspt,even_mons_n1,absorb_tbl)
  local i,j,k,ucoeffs,uidx,uc,err_list,err,m,bad_m,n,c,cm,cp,cfs,sqs,n1,m1,m2;
  printf("\nu = "); lprint(evalf(u));
  ucoeffs,uidx := support(u,X,ms);
  n := nops(X): #bm := [seq(0,i=1..n)]: bp := bm: 
  #for i from 1 to ceil(n/2) do bm[i]:=1:od:for i from ceil(n/2)+1 to n do bp[i]:=1:od:
  #card_absorbed := []; # counts how many monomials to absorbe for each even monomial with support in even_mons
  # printf("\neven_mons = "); lprint(even_mons);
  err_list := [];
  for j from 1 to nops(even_mons) do
    m := even_mons[j];
    #lprint(m);
    n1 := even_mons_n1[j];
    #lprint(n1);
    err := e;
    for i from 1 to nops(ucoeffs) do
      # printf("0 "); lprint(bad_m); lprint(err);
      # printf("1 "): lprint(bad_m); lprint(err);
      uc := ucoeffs[i]; bad_m := mspt[uidx[i]]:
      if m = bad_m then err := err + uc: else if has(bad_m,n1) then err := err - 1/2*abs(uc): fi:fi:
    od:
    err_list := [op(err_list),err];
  od:
  printf("\nerr list = "); lprint(evalf(err_list));
  cfs := [op(err_list)]; sqs := [seq(mul(X[i]^(m[i]/2),i=1..n), m in even_mons)];
  for i from 1 to nops(ucoeffs) do
    uc := ucoeffs[i]; bad_m := mspt[uidx[i]]:
    if not_even(bad_m) then
    #c,cm,cp := decomp_mon(bad_m,n); 
    #printf("\n"); lprint(bad_m);
    #lprint(c+cm); lprint(c+cp);
    #m1 := c+cm: m2 := c + cp:
    k := Search(bad_m,mspt):m1:=mspt[absorb_tbl[k]]:m2:=bad_m - m1;
    printf("\n"); lprint(bad_m);
    lprint(m1); lprint(m2);
    cfs := [op(cfs),abs(uc)/2]; 
    sqs := [op(sqs),mul(X[i]^m1[i],i=1..n)+sign(uc)*mul(X[i]^m2[i],i=1..n)]: fi:
  od:
  return cfs,sqs:
end;

relaxordermin := proc(g)
 return ceil(max(seq(degree(gi),gi in g))/2):
end:

locmatsizes := proc(g,k,n)
  return [seq(binomial(n + k - ceil(degree(gi)/2),n),gi in g)]:
end;

perturbate := proc(p,prec)
  local X,d,k,n,mspt,card_nk,msptk,ms,msk;
  X:=[op(indets(p))]:
  d := degree(p): k := ceil(d/2);
  n := nops(X): 
  mspt := monspt(n,d):
  card_nk := binomial(n+k,n);
  msptk := mspt[1..card_nk];
  mspt,msptk := truncate_withNP(p,X,mspt,msptk):
  ms := mons(n,d,X,mspt): msk := mons(n,k,X,msptk):
  return p + 2^(-prec)*dense_perturbation(msk):
end:

multivsos1:=proc(f,prec::integer := 10, precSVD::integer := 10, precSDP::integer := 200, epsStar::integer := 30, precRound::integer:=30, dataprec::integer:=10, gmp::boolean:=false, algo::integer:=1, glist::list:=[],relaxorder::integer:=0)
  local p,d,mspt,ms,rc,ridx,ng,gc,gidx,gic,giidx,S,s,c, q,n,k,t,e,r,l,a,s1,s2,u,v,i,j,sqs,cfs,sos,rfloat, eigs, eigsg, eigsgi, soslist,soslistg, soslistgi, sumsos,cnd,maxq,obj_plus_r0, card_nk, even_mons,err_list,err,msptk,msk,absorb_tbl,rmin,nsdponeg,idx,oneg,idxi,g,X,ti,tf,gi,cg;
  ti := time[real]():
  c := max(map(_c -> abs(_c),coeffs(expand(f))));
  p := 1/c*f;
  X := [op(indets(p))]:
  g := []; ng := nops(glist);
  for i from 1 to ng do
    gi := expand(glist[i]);
    cg := max(map(_c -> abs(_c),coeffs(expand(gi)))):
    g := [op(g),gi/cg]
  od:
  #g := expand(glist);

  if lcoeff(p) < 0 and ng = 0 then lprint(p): error "There is no decomposition into sum of squares for this polynomial"; fi;

  # Let d = 2 k be the degree; fail if this is odd (can't be pos def) 
  # If it's a constant (i.e. the polynomial is a constant multiple of 
  # a perfect square) then return almost immediately.                 
  if ng = 0 then d := degree(p): else rmin := relaxordermin([p,op(g)]): d := 2*max(relaxorder,rmin): fi:
  k := ceil(d/2);
  if (2 * k <> d and ng = 0) then lprint(p): error "There is no decomposition into sum of squares for this polynomial";  fi;
  if (d = 0 and ng = 0) then lprint(p, " * (", s, ")^2"); fi;
  n := nops(X): 
  mspt := monspt(n,d):
  card_nk := binomial(n+k,n);
  msptk := mspt[1..card_nk];
  if ng = 0 then mspt,msptk := truncate_withNP(p,X,mspt,msptk):fi: # if no constraints then compute NP
  ms := mons(n,d,X,mspt): msk := mons(n,k,X,msptk):


  # let r = p - e * t be a safe perturbation 
  if algo = 1 then
  t := dense_perturbation(msk);
  # printf("\nperturber = "); lprint(t):
  even_mons := 2*msptk;
  # even_mons_n1 := even_n1(n,even_mons);
  e := 1/2^prec: 
  else e := 0: t := 0: 
  fi:

  r := expand(p - e*t);
  rc,ridx := support(r,X,ms):
  gc := []:
  gidx := []:

  for i from 1 to ng do 
    gic,giidx := support(g[i],X,ms):
    gc:=[op(gc),gic]; gidx:=[op(gidx),giidx];
  od:
  (eigs,soslist,obj_plus_r0,absorb_tbl,nsdponeg) := sos2sdp(r,X,k,mspt,msk,rc,ridx,precSVD,precSDP,epsStar,precRound,dataprec,gmp,algo,g,gc,gidx);

  idx := 0: sumsos := obj_plus_r0: oneg := [1,op(g)];
  for i from 1 to ng+1 do
    idxi := idx + nsdponeg[i]:
    sumsos := sumsos + oneg[i]*sum(eigs[j]*soslist[j]^2,j=idx+1..idxi);
    idx := idxi:
  od:
#  sumsos := obj_plus_r0 + sum(eigs[j]*soslist[j]^2,j=1..nops(soslist));
#  for i from 1 to ng do 
#    eigsgi := eigsg[i]: soslistgi := soslistg[i]:
#    sumsos := sumsos + sum(eigsgi[j]*soslistgi[j]^2,j=1..nops(soslistgi))*g[i];
#  od:

  #lprint(evalf(expand(sumsos + obj_plus_r0)));

  if algo = 1  then
    ## if ng > 0 then
    ## t := dense_perturbation(msk);
    ## printf("\nperturber = "); lprint(t):
    ## even_mons := 2*msptk;e := 1/2^prec: r := expand(p - e*t);
    ## fi:
    u := expand(r - sumsos+obj_plus_r0);
    lprint(soslist); lprint(u);
    cfs,sqs := absorber(u,X,e,even_mons,ms,mspt,absorb_tbl);
    nsdponeg[1]:=nsdponeg[1]+nops(sqs):
  else
    u:=0;
    cfs := []: sqs := []:
  fi:
  err := expand(u+e*t - add(cfs[i]*sqs[i]^2,i=1..nops(sqs)));
  if algo = 1 then printf("\nerr = "); lprint(evalf(err)); fi:
  cfs := [op(cfs), op(eigs)]: sqs := [op(sqs),op(soslist)]:
  sos := [];
  for i from 1 to nops(sqs) do
    sos := [op(sos), c*cfs[i],sqs[i]]:
  od;
  
  if soscheck1(f,sos,nsdponeg,g) = 0 then printf("\n Certified SOS Decomposition\n"): fi:
  tf := time[real]()-ti; 
  printf("bitsize= %d\n",BitSizePolSeq(sos,X));
  printf("time= %esecs\n",tf);
  return sos,nsdponeg:
end;

soscheck1:=proc(f, sos, nsdponeg, g::list:=[])
  local s,i,j,idx,idxi,oneg;
  oneg := [1,op(g)]:
  s := 0; idx:=0;
for i from 1 to nops(g)+1 do
  idxi:=idx+nsdponeg[i]:
  for j from idx+1 to idxi do 
    if sos[2*j-1] < 0 then lprint(evalf(sos[2*j-1])): error "Negative number => invalid sum of squares decomposition"; 
    else  s := s + oneg[i]*sos[2*j-1]*sos[2*j]^2 fi: 
  od:
  idx:=idxi:
od:
  if not expand(f - s) = 0 then lprint(evalf(expand(f - s))); error "Inexact sum of squares decomposition"; else return 0: fi:
end;

BitSizeSos := proc (sos,X)
   return add(BitSizePol(p[1], X) + BitSizePolQuadr(p[2], X), p in sos);
end;

BitSizePolQuadr := proc(q,X)
  return BitRat(q[1]) + BitSizePol(q[2], X) + BitRat(q[3]);
end;

BitSizePolSeq3 := proc(listpol,X)
  return add(BitSizePol(p[1], X) + BitSizePolSeq(p[2], X), p in listpol);
end;

BitSizePolSeq := proc(listpol,X)
  return add(BitSizePol(p, X), p in listpol);
end;

BitSizePol:=proc(p, X)
  local res;
  res := [coeffs(expand(p),X)];
  return BitSizeSeq(res);
end;

BitSizeSeq:=proc(l)
  return add(BitRat(c), c in l);
end;

BitRat := proc(r)
  local n, d, res,rs;
  if type(r,rational) then rs :=r : else rs:=r^2: fi:
  if rs = 0 then return 1; fi;
  (n, d) := (abs(numer(rs)), abs(denom(rs)));
  if d = 1 then res :=  ilog2(n) + 1 else res := ilog2(n) + ilog2(d) + 2 fi;
  return res;
end;

benchRAGLib := proc(f,g::list:=[])
  local sys,i,ti,tf,sol;
  ti := time[real]():
  sys := [expand(f) < 0];
  for i from 1 to nops(g) do
    sys := [op(sys), expand(-g[i]) < 0]
  od:
  sol := HasRealSolutions(sys);
  tf := time[real]()-ti; 
  printf("time= %esecs\n",tf);
  return sol;
end;

benchSamplePoints := proc(f,g::list:=[])
  local sys,vars,R,P,i,ti,tf;
  ti := time[real]():
  sys := [f < 0];
  for i from 1 to nops(g) do
    sys := [op(sys), g[i] >= 0]
  od:
  vars:=[op(indets(f))]: R:=PolynomialRing(vars): P:=SamplePoints(sys,R);
  tf := time[real]()-ti; 
  printf("time= %esecs\n",tf);
  return P;
end;
