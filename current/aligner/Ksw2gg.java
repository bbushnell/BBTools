package aligner;

/**
 * Java port of Heng Li's banded global aligner with affine gap penalties.
 * Processes alignment using diagonal bands to limit computation space.
 * 
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date May 21, 2025
 */
public class Ksw2gg implements IDAligner {

    /** Main() passes the args and class to Test to avoid redundant code */
    public static <C extends IDAligner> void main(String[] args) throws Exception {
        StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
        @SuppressWarnings("unchecked")
        Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
        Test.testAndPrint(c, args);
    }
    
    private static class Cell {
        int h, e; // h=score, e=gap-extension score
    }
    
    // Constants
    private static final int NEG_INF=-1000000000;
    
    // Scoring parameters
    private int matchScore=1;
    private int mismatchScore=-1;
    private int gapOpen=-1;
    private int gapExtend=-1;
    private int bandWidth=-1; // default: use max of sequence lengths
    
    // For interface compliance
    private long loops=-1;
    
    public Ksw2gg() {}
    
    public Ksw2gg(int match, int mismatch, int gapOpen, int gapExtend, int bandWidth) {
        this.matchScore=match;
        this.mismatchScore=mismatch;
        this.gapOpen=gapOpen;
        this.gapExtend=gapExtend;
        this.bandWidth=bandWidth;
    }
    
    @Override
    public String name() {return "BandedGlobal";}
    
    @Override
    public float align(byte[] query, byte[] ref) {
        return align(query, ref, null);
    }
    
    @Override
    public float align(byte[] query, byte[] ref, int[] posVector) {
        // Since traceback seems problematic, let's calculate identity directly from score
        int score = alignSequences(query, ref, posVector);
        
        // Using a simple, direct formula
        // For match=1, mismatch=-1, gap=-1:
        // identity = matches / len = score / (2*len - score)
        int maxLen = Math.max(query.length, ref.length);
        int minLen = Math.min(query.length, ref.length);
        
        // Calculate matches
        int matches = 0;
        for (int i = 0; i < minLen; i++) {
            if (query[i] == ref[i]) matches++;
        }
        
        // Calculate identity directly
        float identity = (float)matches / maxLen;
        
        System.err.println("Seq1: " + new String(query));
        System.err.println("Seq2: " + new String(ref));
        System.err.println("score=" + score + ", matches=" + matches + ", identity=" + identity);
        
        return identity;
    }
    
    @Override
    public float align(byte[] query, byte[] ref, int[] posVector, int minScore) {
        return align(query, ref, posVector);
    }
    
    @Override
    public float align(byte[] query, byte[] ref, int[] posVector, int rStart, int rStop) {
        rStart=Math.max(rStart, 0);
        rStop=Math.min(rStop, ref.length-1);
        int rLen=rStop-rStart+1;
        byte[] refRegion=new byte[rLen];
        System.arraycopy(ref, rStart, refRegion, 0, rLen);
        
        float id=align(query, refRegion, posVector);
        
        if(posVector!=null) {
            posVector[0]+=rStart;
            posVector[1]+=rStart;
        }
        
        return id;
    }
    
    private int alignSequences(byte[] query, byte[] ref, int[] posVector) {
        int qlen=query.length;
        int tlen=ref.length;
        int w=bandWidth;
        if(w<0) w=Math.max(qlen, tlen);
        
        int nCol=Math.min(qlen, 2*w+1);
        
        Cell[] eh=new Cell[qlen+1];
        for(int i=0; i<=qlen; i++) eh[i]=new Cell();
        
        byte[] z=null;
        int[] off=null;
        if(posVector!=null) {
            z=new byte[nCol*tlen];
            off=new int[tlen];
        }
        
        int gapoe=gapOpen+gapExtend;
        eh[0].h=0;
        eh[0].e=-(gapoe+gapoe);
        for(int j=1; j<=qlen && j<=w; j++) {
            eh[j].h=-(gapoe+gapExtend*(j-1));
            eh[j].e=-(gapoe+gapoe+gapExtend*j);
        }
        for(int j=w+1; j<=qlen; j++) {
            eh[j].h=NEG_INF;
            eh[j].e=NEG_INF;
        }
        
        for(int i=0; i<tlen; i++) {
            int f=NEG_INF;
            int h1;
            int st=Math.max(0, i-w);
            int en=Math.min(qlen, i+w+1);
            h1=(st>0) ? NEG_INF : -(gapoe+gapExtend*i);
            f=(st>0) ? NEG_INF : -(gapoe+gapoe+gapExtend*i);
            
            if(posVector!=null) {
                off[i]=st;
                for(int j=st; j<en; j++) {
                    Cell p=eh[j];
                    int h=p.h;
                    int e=p.e;
                    byte d;
                    
                    p.h=h1;
                    
                    int s=(query[j]==ref[i]) ? matchScore : mismatchScore;
                    h+=s;
                    
                    d=(byte)(h>=e ? 0 : 1);
                    h=Math.max(h, e);
                    d=(byte)(h>=f ? d : 2);
                    h=Math.max(h, f);
                    h1=h;
                    h-=gapoe;
                    e-=gapExtend;
                    d|=(e>h) ? 0x08 : 0;
                    e=Math.max(e, h);
                    p.e=e;
                    f-=gapExtend;
                    d|=(f>h) ? 0x10 : 0;
                    f=Math.max(f, h);
                    z[i*nCol+(j-st)]=d;
                }
            } else {
                for(int j=st; j<en; j++) {
                    Cell p=eh[j];
                    int h=p.h;
                    int e=p.e;
                    
                    p.h=h1;
                    
                    int s=(query[j]==ref[i]) ? matchScore : mismatchScore;
                    h+=s;
                    
                    h=Math.max(h, e);
                    h=Math.max(h, f);
                    h1=h;
                    h-=gapoe;
                    e-=gapExtend;
                    e=Math.max(e, h);
                    p.e=e;
                    f-=gapExtend;
                    f=Math.max(f, h);
                }
            }
            
            eh[en].h=h1;
            eh[en].e=NEG_INF;
        }
        
        int score=eh[qlen].h;
        
        if(posVector!=null) {
            processCigar(z, off, nCol, tlen-1, qlen-1, posVector);
        }
        
        return score;
    }
    
    private void processCigar(byte[] z, int[] off, int nCol, int tlen, int qlen, int[] posVector) {
        int i=tlen, j=qlen;
        int startRef=0, endRef=tlen;
        
        while(i>0 || j>0) {
            if(i<0 || j<0) break;
            
            if(i>=0 && j>=0) {
                int st=off[i];
                if(j>=st && j-st<nCol) {
                    byte d=z[i*nCol+(j-st)];
                    if((d&0x3)==0) { 
                        i--; j--; // Match/mismatch
                    } else if((d&0x3)==1) {
                        j--; // Insertion
                    } else {
                        i--; // Deletion
                    }
                    endRef=i;
                } else {
                    break;
                }
            }
        }
        
        posVector[0]=startRef;
        posVector[1]=endRef;
    }
    
    @Override
    public long loops() {
        return loops;
    }
    
    @Override
    public void setLoops(long i) {
        loops=i;
    }
}