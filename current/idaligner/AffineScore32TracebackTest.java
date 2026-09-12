package idaligner;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/** Independent exhaustive edit-path oracle plus path replay and dense-score comparisons. */
public final class AffineScore32TracebackTest {
 public static void main(String[] args) {
  boolean enabled=false; assert(enabled=true);
  if(!enabled){throw new IllegalStateException("Assertions required");}
  int[][] profiles={{100,-127,0,-395,-39,-472,-33},{5,-4,-1,-7,-2,-9,-1},{3,0,3,0,0,0,0},{7,-2,9,-1,-9,-8,0}};
  List<byte[]> small=new ArrayList<byte[]>();
  words(small,"",3);
  Random random=new Random(20260912570L);
  int exhaustive=0, randomized=0, windows=0;
  System.out.println("profile\texhaustive_paths_cases\trandom_dense_cases\twindow_cases\tpath_rescore");
  for(int k=0;k<profiles.length;k++){
   int[] p=profiles[k];
   AffineScore32Traceback tracer=tracer(p); AffineScore32 dense=dense(p);
   for(byte[] q:small){for(byte[] r:small){
    AffineScore32Traceback.Result result=tracer.align(q,r);
    check(q,r,0,r.length,result,p);
    assert(result.score==exhaustive(q,r,p)) : "Independent edit-path optimum mismatch";
    assert(result.score==dense.score(q,r)) : "Score-only API mismatch";
    exhaustive++;
   }}
   for(int t=0;t<160;t++){
    byte[] q=bases(random,random.nextInt(33)), r=bases(random,1+random.nextInt(48));
    byte[] qCopy=q.clone(), rCopy=r.clone();
    AffineScore32Traceback.Result result=tracer.align(q,r), again=tracer.align(q,r);
    check(q,r,0,r.length,result,p);
    assert(result.score==dense.score(q,r)) : "Random dense optimum mismatch";
    assert(result.referenceStart==again.referenceStart && result.referenceEndExclusive==again.referenceEndExclusive && Arrays.equals(result.operations(),again.operations())) : "Traceback ties must be deterministic";
    assert(Arrays.equals(q,qCopy) && Arrays.equals(r,rCopy)) : "Scoring must not mutate caller input";
    byte[] padded=new byte[r.length+6]; Arrays.fill(padded,(byte)'N'); System.arraycopy(r,0,padded,3,r.length);
    AffineScore32Traceback.Result window=tracer.align(q,padded,3,r.length+2);
    check(q,padded,3,r.length+3,window,p);
    assert(window.score==result.score && window.referenceStart==result.referenceStart+3 && window.referenceEndExclusive==result.referenceEndExclusive+3 && Arrays.equals(window.operations(),result.operations())) : "Window translation changes alignment";
    byte[] saved=result.operations(); tracer.align(bytes("ACGTACGT"),bytes("TTACGTACGTTT"));
    assert(Arrays.equals(saved,result.operations())) : "Later calls overwrite retained results";
    if(saved.length>0){byte[] external=result.operations(); external[0]='?'; assert(Arrays.equals(saved,result.operations())) : "Result operation array leaked";}
    randomized++; windows++;
   }
   System.out.println(k+"\t1600\t160\t160\tEXACT");
  }
  AffineScore32Traceback allI=new AffineScore32Traceback(1,-100,0,0,0,-100,-100);
  AffineScore32Traceback.Result gap=allI.align(bytes("AA"),bytes("C"));
  assert(gap.score==0 && gap.referenceStart==1 && gap.referenceEndExclusive==1 && Arrays.equals(gap.operations(),bytes("II"))) : "All-insertion endpoint must preserve empty reference span";
  AffineScore32Traceback.Result leftmost=AffineScore32Traceback.msaLike().align(bytes("A"),bytes("AAA"));
  assert(leftmost.score==100 && leftmost.referenceStart==0 && leftmost.referenceEndExclusive==1 && Arrays.equals(leftmost.operations(),bytes("M"))) : "Equal perfect placements must choose the earliest final column";
  AffineScore32Traceback.Result stateTie=new AffineScore32Traceback(3,0,3,0,0,0,0).align(bytes("AA"),bytes("A"));
  assert(stateTie.score==3 && stateTie.referenceStart==0 && stateTie.referenceEndExclusive==1 && Arrays.equals(stateTie.operations(),bytes("IM"))) : "Final M wins the score tie with I; its predecessor is the leading insertion";
  int rejects=0;
  try {AffineScore32Traceback.msaLike().align(new byte[50000],new byte[50000]);} catch(IllegalArgumentException expected){rejects++;}
  try {new AffineScore32Traceback(Integer.MAX_VALUE,-1,0,-1,-1,-1,-1).align(bytes("AA"),bytes("A"));} catch(IllegalArgumentException expected){rejects++;}
  try {AffineScore32Traceback.msaLike().align(bytes("A"),bytes("A"),1,0);} catch(IndexOutOfBoundsException expected){rejects++;}
  assert(rejects==3) : "Reject unsafe trace products, scores and reference windows before allocation";
  AffineScore32Traceback memory=AffineScore32Traceback.msaLike();
  memory.align(new byte[100],new byte[150]);
  assert(memory.retainedTraceBytes()==101*151 && memory.retainedScoreBytes()==6*4*151) : "Traceback storage contract changed";
  assert(exhaustive==6400 && randomized==640 && windows==640) : "Coverage counters must match the declared panel";
  System.out.println("PASS exhaustive="+exhaustive+" random="+randomized+" windows="+windows+" ties=2 rejects="+rejects+" traceBytes="+memory.retainedTraceBytes()+" scoreBytes="+memory.retainedScoreBytes());
 }
 private static void check(byte[] q,byte[] r,int start,int end,AffineScore32Traceback.Result a,int[] p){
  assert(a.referenceStart>=start && a.referenceEndExclusive<=end && a.referenceStart<=a.referenceEndExclusive) : "Reference endpoint outside requested window";
  int qi=0,ri=a.referenceStart; long score=0; byte previous='M';
  for(byte op:a.operations()){
   if(op=='M'){assert(qi<q.length && ri<a.referenceEndExclusive) : "Diagonal exceeds consumed spans"; score+=pair(q[qi++],r[ri++],p);}
   else if(op=='I'){assert(qi<q.length) : "Insertion exceeds query"; qi++; score+=previous=='I'?p[4]:p[3];}
   else {assert(op=='D' && ri<a.referenceEndExclusive) : "Invalid deletion/path code"; ri++; score+=previous=='D'?p[6]:p[5];}
   previous=op;
  }
  assert(qi==q.length && ri==a.referenceEndExclusive && score==a.score) : "Path consumption/rescore differs from reported optimum: "+score+" != "+a.score;
 }
 private static long exhaustive(byte[] q,byte[] r,int[] p){
  if(q.length==0){return 0;}
  if(r.length==0){return (long)p[3]+(q.length-1L)*p[4];}
  long best=Long.MIN_VALUE;
  for(int start=0;start<=r.length;start++){best=Math.max(best,walk(q,r,p,0,start,'M',0));}
  assert(best>Long.MIN_VALUE) : "A complete-query path must exist"; return best;
 }
 private static long walk(byte[] q,byte[] r,int[] p,int qi,int ri,char prev,long score){
  assert(qi<=q.length && ri<=r.length) : "Enumerated path left its finite grid";
  long best=qi==q.length && ri>0 ? score : Long.MIN_VALUE;
  if(qi<q.length){
   best=Math.max(best,walk(q,r,p,qi+1,ri,'I',score+(prev=='I'?p[4]:p[3])));
   if(ri<r.length){best=Math.max(best,walk(q,r,p,qi+1,ri+1,'M',score+pair(q[qi],r[ri],p)));}
  }
  if(ri<r.length && qi>0){best=Math.max(best,walk(q,r,p,qi,ri+1,'D',score+(prev=='D'?p[6]:p[5])));}
  return best;
 }
 private static int pair(byte a,byte b,int[] p){
  char x=Character.toUpperCase((char)a), y=Character.toUpperCase((char)b); if(x=='U'){x='T';} if(y=='U'){y='T';}
  return "ACGT".indexOf(x)<0 || "ACGT".indexOf(y)<0 ? p[2] : x==y ? p[0] : p[1];
 }
 private static void words(List<byte[]> out,String word,int remaining){out.add(bytes(word)); if(remaining>0){words(out,word+"A",remaining-1); words(out,word+"C",remaining-1); words(out,word+"N",remaining-1);}}
 private static byte[] bases(Random r,int n){byte[] b=new byte[n]; for(int i=0;i<n;i++){b[i]=(byte)"ACGTNacu".charAt(r.nextInt(8));} return b;}
 private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
 private static AffineScore32Traceback tracer(int[] p){return new AffineScore32Traceback(p[0],p[1],p[2],p[3],p[4],p[5],p[6]);}
 private static AffineScore32 dense(int[] p){return new AffineScore32(p[0],p[1],p[2],p[3],p[4],p[5],p[6]);}
}
