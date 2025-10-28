package stream;

import java.util.Arrays;

/**
 * Minimal driver to exercise MDWalker and MDWalker2 on synthetic cases.
 * Prints initial/updated longmatch for visual inspection.
 */
public class TestMDWalker {

    private static byte[] lm(String s){return s.getBytes();}

    private static void runCase(String name, String md, String cigar, String longmatch){
        System.out.println("== "+name+" ==");
        byte[] lm0=lm(longmatch);
        byte[] lm1=Arrays.copyOf(lm0, lm0.length);
        byte[] lm2=Arrays.copyOf(lm0, lm0.length);

        try{
            MDWalker w=new MDWalker(md, cigar, lm1, null);
            w.fixMatch(null);
            System.out.println("MDWalker : "+new String(lm1));
        }catch(Throwable t){
            System.out.println("MDWalker : threw "+t.getClass().getSimpleName()+" - "+t.getMessage());
        }

        try{
            MDWalker2 w2=new MDWalker2(md, cigar, lm2, null);
            w2.fixMatch(null);
            System.out.println("MDWalker2: "+new String(lm2));
        }catch(Throwable t){
            System.out.println("MDWalker2: threw "+t.getClass().getSimpleName()+" - "+t.getMessage());
        }

        System.out.println("Start    : "+longmatch);
        System.out.println();
    }

    public static void main(String[] args){
        // 1) Simple substitution after 3 matches: MD:Z:3A2
        runCase("simple-sub", "MD:Z:3A2", "6M", "mmmmmm");

        // 2) Insertion before substitution should be skipped by match counter
        // CIGAR: 2M1I3M, MD: 2A3
        runCase("ins-before-sub", "MD:Z:2A3", "2M1I3M", "mmImmm");

        // 3) Deletion before substitution: 2M2D1M1T2M → MD: 2^CC1T2
        runCase("del-before-sub", "MD:Z:2^CC1T2", "2M2D1M1X2M", "mmDDmSmm");

        // 4) Leading clipping should be skipped safely
        runCase("leading-clip", "MD:Z:1A1", "1S1M1X1M", "C m S m".replace(" ",""));

        // 5) Mixed I/D around matches
        runCase("mixed-indels", "MD:Z:1A1^G2", "1M1I1M1D2M", "mImDmm");
    }
}

