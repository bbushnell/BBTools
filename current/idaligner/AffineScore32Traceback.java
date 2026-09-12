package idaligner;

import java.util.Arrays;

/** Experimental traceback for AffineScore32's affine model, not MSA scoring.
 * Full query, free reference ends. Six reusable score rows and one predecessor
 * byte per grid cell; no SIMD/pruning or mapper integration is provided here.
 * Each worker owns an instance. Results own their operation arrays.
 */
public final class AffineScore32Traceback {
 public AffineScore32Traceback(int match, int substitution, int ambiguity,
   int insertionOpen, int insertionExtend, int deletionOpen, int deletionExtend) {
  model=new AffineScore32(match,substitution,ambiguity,insertionOpen,insertionExtend,deletionOpen,deletionExtend);
 }
 public static AffineScore32Traceback msaLike() {
  return new AffineScore32Traceback(100,-127,0,-395,-39,-472,-33);
 }
 public Result align(byte[] query, byte[] reference) {
  if(query==null || reference==null){throw new NullPointerException("Query and reference are required");}
  if(reference.length==0){
   int score=model.scoreEmptyReference(query.length);
   byte[] ops=new byte[query.length]; Arrays.fill(ops,(byte)'I');
   return new Result(score,0,0,ops);
  }
  return align(query,reference,0,reference.length-1);
 }
 /** Inclusive nonempty input window; result coordinates are absolute, half-open.
  * An all-insertion path can consume no reference: start==endExclusive.
  * Ties prefer the earliest final reference column, then states M,I,D;
  * predecessor ties also prefer M,I,D. This is not an MSA tie convention.
  */
 public Result align(byte[] query, byte[] reference, int start, int end) {
  if(query==null || reference==null){throw new NullPointerException("Query and reference are required");}
  if(start<0 || end<start || end>=reference.length){throw new IndexOutOfBoundsException("Invalid inclusive reference window");}
  int q=query.length, r=end-start+1;
  AffineScore32.validateScoreRange(q,r,model.match,model.substitution,model.ambiguity,
    model.insertionOpen,model.insertionExtend,model.deletionOpen,model.deletionExtend);
  if(q==0){return new Result(0,start,start,new byte[0]);}
  long cells=(q+1L)*(r+1L);
  if(cells>Integer.MAX_VALUE-8L){throw new IllegalArgumentException("Trace grid exceeds one byte-array capacity: "+cells);}
  ensureCapacity(r+1,(int)cells);
  int[] pm=rows[0], pi=rows[1], pd=rows[2], cm=rows[3], ci=rows[4], cd=rows[5];
  Arrays.fill(pm,0,r+1,0); Arrays.fill(pi,0,r+1,BAD); Arrays.fill(pd,0,r+1,0);
  int stride=r+1;
  for(int i=1;i<=q;i++){
   cm[0]=cd[0]=BAD;
   ci[0]=(i==1 ? model.insertionOpen : add(pi[0],model.insertionExtend));
   trace[i*stride]=(byte)((i==1 ? M : I)<<2);
   for(int j=1;j<=r;j++){
    int mp=winner(pm[j-1],pi[j-1],pd[j-1]);
    cm[j]=add(value(mp,pm[j-1],pi[j-1],pd[j-1]),model.baseScore(query[i-1],reference[start+j-1]));
    int im=add(pm[j],model.insertionOpen), ii=add(pi[j],model.insertionExtend), id=add(pd[j],model.insertionOpen);
    int ip=winner(im,ii,id); ci[j]=value(ip,im,ii,id);
    int dm=add(cm[j-1],model.deletionOpen), di=add(ci[j-1],model.deletionOpen), dd=add(cd[j-1],model.deletionExtend);
    int dp=winner(dm,di,dd); cd[j]=value(dp,dm,di,dd);
    trace[i*stride+j]=(byte)(mp | (ip<<2) | (dp<<4));
   }
   int[] swap=pm; pm=cm; cm=swap; swap=pi; pi=ci; ci=swap; swap=pd; pd=cd; cd=swap;
  }
  int best=BAD, column=-1, state=-1;
  // AffineScore32 excludes final column0 when the reference window is nonempty.
  for(int j=1;j<=r;j++){
   int s=winner(pm[j],pi[j],pd[j]), score=value(s,pm[j],pi[j],pd[j]);
   if(score>best){best=score; column=j; state=s;}
  }
  assert(column>=1 && state>=M && state<=D && best>BAD) : "Numeric bounds must leave a reachable endpoint";
  return decode(q,column,state,stride,start,best);
 }
 private Result decode(int row, int column, int state, int stride, int start, int score) {
  assert(row>0 && column>0) : "Nonempty-query traceback starts at a reachable final-row column";
  int finish=start+column, count=0;
  byte[] reverse=new byte[row+column];
  while(row>0){
   assert(column>=0 && state>=M && state<=D) : "Predecessor must stay inside the affine grid";
   int previous=(trace[row*stride+column] >>> (state*2))&3;
   if(state==M){assert(column>0) : "Diagonal consumes a reference base"; reverse[count++]='M'; row--; column--;}
   else if(state==I){reverse[count++]='I'; row--;}
   else {assert(column>0) : "Deletion consumes a reference base"; reverse[count++]='D'; column--;}
   state=previous;
  }
  byte[] operations=new byte[count];
  for(int i=0;i<count;i++){operations[i]=reverse[count-i-1];}
  return new Result(score,start+column,finish,operations);
 }
 private void ensureCapacity(int columns, int cells) {
  assert(columns>1 && cells>=columns) : "Trace dimensions must include the boundary row/column";
  if(rows[0].length<columns){for(int i=0;i<rows.length;i++){rows[i]=new int[columns];}}
  if(trace.length<cells){trace=new byte[cells];}
 }
 private static int winner(int m,int i,int d){return m>=i && m>=d ? M : i>=d ? I : D;}
 private static int value(int state,int m,int i,int d){assert(state>=M && state<=D) : "Invalid affine state"; return state==M ? m : state==I ? i : d;}
 private static int add(int score,int delta){return score<=BAD ? BAD : score+delta;}
 /** Retained primitive payload only; excludes object headers and per-call operations/scratch. */
 public long retainedTraceBytes(){return trace.length;}
 /** Retained primitive score-row payload, not total heap or process memory. */
 public long retainedScoreBytes(){return 6L*Integer.BYTES*rows[0].length;}

 public static final class Result {
  private Result(int score_,int start_,int end_,byte[] ops){
   assert(start_>=0 && end_>=start_) : "Half-open aligned reference span cannot be reversed";
   score=score_; referenceStart=start_; referenceEndExclusive=end_; operations=ops;
  }
  /** M consumes one query/reference pair (including mismatch or ambiguity); I/D consume query/reference respectively. */
  public byte[] operations(){return operations.clone();}
  public final int score, referenceStart, referenceEndExclusive;
  private final byte[] operations;
 }
 private final AffineScore32 model;
 private final int[][] rows={new int[0],new int[0],new int[0],new int[0],new int[0],new int[0]};
 private byte[] trace=new byte[0];
 private static final int M=0,I=1,D=2,BAD=AffineScore32.BAD;
}
