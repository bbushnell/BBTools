package cardinality;

import shared.Tools;

/**
 * Gemini's rewrite of MantissaCompare. Single-bucket 16-bit simulation.
 * State 0 = Established/Small Gap. State Cap = Outlier/Big Gap.
 */
public class MantissaCompare2 {

    static final int MODE_MANTISSA=0, MODE_ANDTISSA=1, MODE_NLZ2=2, MODE_HISTORY=3, MODE_LUCK=4;
    static final String[] MODE_NAMES={"Mantissa", "Andtissa", "NLZ2", "History", "Luck"};

    public static void main(String[] args){
        int inner=32768, outer=131072, maxTier=11, mode=MODE_MANTISSA, bits=2;

        for(String arg : args){
            final String[] ab=arg.split("=");
            final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
            if(a.equals("inner")){inner=Integer.parseInt(b);}
            else if(a.equals("outer")){outer=Integer.parseInt(b);}
            else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
            else if(a.equals("bits")){bits=Integer.parseInt(b);}
            else if(a.equals("mode")){
                if(b.equals("mantissa")){mode=MODE_MANTISSA;}
                else if(b.equals("andtissa")){mode=MODE_ANDTISSA;}
                else if(b.equals("nlz2")){mode=MODE_NLZ2;}
                else if(b.equals("history")){mode=MODE_HISTORY;}
                else if(b.equals("luck")){mode=MODE_LUCK;}
            }
        }

        final int numStates=1<<bits;
        final int cap=numStates-1;
        System.err.println("MantissaCompare2: mode="+MODE_NAMES[mode]+" bits="+bits
            +" states="+numStates+" inner="+inner+" outer="+outer+" maxtier="+maxTier);
        final long t0=System.nanoTime();

        final long[][] count=new long[64][numStates];
        final double[][] sum=new double[64][numStates];

        for(int trial=0; trial<outer; trial++){
            long seed=Tools.hash64shift((long)trial*1234567891L+54321);
            char stored=0;
            int maxNlz=-1, extra=0, luckSecond=-1;

            for(int card=1; card<=inner; card++){
                final long key=Tools.hash64shift(seed+card);
                final int nlz=Long.numberOfLeadingZeros(key);

                if(mode==MODE_MANTISSA || mode==MODE_ANDTISSA || mode==MODE_NLZ2){
                    int val=0;
                    if(mode==MODE_MANTISSA){
                        final int shift=63-nlz-bits;
                        val=(shift<0) ? 0 : (int)((~(key>>>shift))&cap);
                    }else if(mode==MODE_ANDTISSA){
                        final long anded=key&(key<<2);
                        final int shift=63-nlz-bits;
                        val=(shift<0) ? 0 : (int)((anded>>>shift)&cap);
                    }else{
                        final int consumed=nlz+1;
                        val=(consumed>=64) ? cap : Math.min(cap, Long.numberOfLeadingZeros(key<<consumed));
                    }
                    if(nlz>maxNlz){ maxNlz=nlz; extra=val; }
                    else if(nlz==maxNlz && val>extra){ extra=val; }

                }else if(mode==MODE_HISTORY){
                    final int newNlzStored=Math.min(nlz+1, 63);
                    final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;
                    if(newNlzStored>oldNlzStored){
                        final int k=newNlzStored-oldNlzStored;
                        final int oldHist=(stored>0) ? (stored&cap) : 0;
                        final int carry=(stored>0) ? (1<<bits) : 0;
                        extra=((oldHist|carry)>>k)&cap;
                        maxNlz=newNlzStored-1;
                    }else if(newNlzStored<oldNlzStored){
                        final int diff=oldNlzStored-newNlzStored;
                        if(diff>=1 && diff<=bits){ extra|=(1<<(bits-diff)); }
                    }
                }else if(mode==MODE_LUCK){
                    if(nlz>maxNlz){ luckSecond=maxNlz; maxNlz=nlz; }
                    else if(nlz>luckSecond && nlz<maxNlz){ luckSecond=nlz; }
                    extra=(maxNlz>=0 && luckSecond>=0) ? Math.min(cap, Math.max(0, maxNlz-luckSecond-1)) : cap;
                }

                if(maxNlz>=0){
                    stored=(char)(((maxNlz+1)<<10)|extra);
                    if(maxNlz<=maxTier && extra<numStates){ count[maxNlz][extra]++; sum[maxNlz][extra]+=card; }
                }
            }
        }

        final double elapsed=(System.nanoTime()-t0)*1e-9;
        System.err.println(String.format("Elapsed: %.1fs", elapsed));

        // Report per-tier
        StringBuilder header=new StringBuilder("Tier\tTotal");
        for(int s=0; s<numStates; s++){header.append("\tP("+s+")\tCF("+s+")");}
        System.out.println(header);

        for(int t=0; t<=maxTier; t++){
            long tierTotal=0; double tierSum=0;
            for(int s=0; s<numStates; s++){ tierTotal+=count[t][s]; tierSum+=sum[t][s]; }
            if(tierTotal<100) continue;
            double tierAvg=tierSum/tierTotal;
            StringBuilder sb=new StringBuilder();
            sb.append(t).append('\t').append(tierTotal);
            for(int s=0; s<numStates; s++){
                sb.append('\t').append(String.format("%.4f", count[t][s]/(double)tierTotal));
                if(count[t][s]>10){
                    sb.append('\t').append(String.format("%+.4f", Math.log(sum[t][s]/count[t][s]/tierAvg)/Math.log(2)));
                }else{sb.append("\tN/A");}
            }
            System.out.println(sb);
        }

        // Weighted summary
        long gTotal=0; double[] gCF=new double[numStates];
        for(int t=3; t<=maxTier; t++){
            long tt=0; double ts=0;
            for(int s=0; s<numStates; s++){ tt+=count[t][s]; ts+=sum[t][s]; }
            if(tt<100) continue;
            double ta=ts/tt;
            for(int s=0; s<numStates; s++){
                if(count[t][s]>10){ gCF[s]+=Math.log(sum[t][s]/count[t][s]/ta)/Math.log(2)*count[t][s]; }
                gTotal+=count[t][s];
            }
        }
        System.out.println("\n# "+MODE_NAMES[mode]+" "+bits+"-bit weighted CFs (tiers 3-"+maxTier+"):");
        StringBuilder cfLine=new StringBuilder("# CF=[");
        for(int s=0; s<numStates; s++){
            if(s>0) cfLine.append(", ");
            cfLine.append(String.format("%+.6f", gCF[s]/(gTotal>0?gTotal:1)));
        }
        cfLine.append("]");
        System.out.println(cfLine);
    }
}
