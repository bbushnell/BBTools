package cardinality;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import shared.Timer;
import shared.Tools;

/**
 * Evaluates cardinality estimators on a REAL (anonymized) stream: a text file
 * of A48-encoded signed longs, one per line — the format 45's re-encoder
 * emits.  The values are already anonymized hashcodes; estimators rehash
 * internally, which is harmless (hash of hash is uniform).  Ground truth is
 * tracked exactly with an open-addressing long set, and estimates are
 * reported at every doubling of distinct plus at end-of-stream.
 *
 * Lanes: AVLL (equal bytes), FLL53 MLE, FLL53 blend (snapshot net), and
 * FLL53 chron (two-clock trajectory fusion) when a chron net is given.
 *
 * Usage: java cardinality.HashStreamEval in=<file> [orgs=384]
 *        [table=fll53mle_384o.tsv] [net=netwin53.txt] [chron=netchron53c.txt]
 *
 * @author Amber
 * @date July 2026
 */
public class HashStreamEval {

	public static void main(String[] args) throws Exception{
		final Timer t=new Timer();
		String in=null, table=null, net="netwin53.txt", chron=null;
		int orgs=384;
		boolean guard=false;
		for(String arg : args){
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			final String b=(split.length>1) ? split[1] : null;
			if(a.equals("in")){in=b;}
			else if(a.equals("orgs")){orgs=Integer.parseInt(b);}
			else if(a.equals("table")){table=b;}
			else if(a.equals("net")){net=b;}
			else if(a.equals("chron")){chron=b;}
			else if(a.equals("guard")){guard=Parse.parseBoolean(b);}
			else if(a.equals("verbose")){/*quiet*/}
			else{throw new IllegalArgumentException("Unknown arg: "+arg);}
		}
		assert(in!=null) : "in=<file> is required";
		if(table==null){table="fll53mle_"+orgs+"o.tsv";}
		CardinalityTracker.clampToAdded=false;

		final int B=orgs*5;
		final int avllRegs=(int)(orgs*8L/3*11/8);   // equal bytes convention
		final ArithmeticVariableLogLog avll=new ArithmeticVariableLogLog(avllRegs, 31, -1, 0);
		final Fll53 fll=new Fll53(B, 31, -1, 0);
		final Fll53MLE mle=new Fll53MLE(table);
		final ComplexityFarm.BigNet win=ComplexityFarm.BigNet.load(net);
		final ComplexityFarm.BigNet fuse=(chron!=null) ? ComplexityFarm.BigNet.load(chron) : null;
		// v5a envelope guard (guard=t): fall back chron->blend when fusion
		// features leave the training envelope (chron+".env", written at train
		// time).  Off by default: old behavior exactly.
		final double[][] env=(guard && chron!=null) ? Fll53Chron.loadEnv(chron+".env") : null;
		if(guard && chron!=null && env==null){
			System.err.println("WARNING: guard=t but "+chron+".env not found; guard inactive");
		}
		final Fll53Chron.Traj tr=new Fll53Chron.Traj(win, mle, 4L*B);
		System.err.println("orgs="+orgs+" bytes="+(orgs*8/3)+" avllRegs="+avllRegs
			+" table="+table+" net="+net+" chron="+chron+" guard="+(env!=null));

		// Exact truth: open-addressing long set (0 sentinel handled separately).
		long[] set=new long[1<<20];
		int setSize=0;
		boolean seenZero=false;

		long adds=0, distinct=0, nextCp=1;
		System.out.println("#adds\tdistinct\tAVLL\tMLE\tblend\tchron\teA%\teM%\teB%\teC%");

		final ByteFile bf=ByteFile.makeByteFile(
			FileFormat.testInput(in, FileFormat.TXT, null, true, true));
		final LineParser1 lp=new LineParser1('\t');
		byte[] line;
		while((line=bf.nextLine())!=null){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			final long v=lp.parseLongA48(0);
			// truth
			if(v==0){
				if(!seenZero){seenZero=true; distinct++;}
			}else{
				int mask=set.length-1;
				int pos=(int)(Tools.hash64shift(v))&mask;
				boolean found=false;
				while(set[pos]!=0){
					if(set[pos]==v){found=true; break;}
					pos=(pos+1)&mask;
				}
				if(!found){
					set[pos]=v; setSize++; distinct++;
					if(setSize*2>set.length){   // grow at 50% load
						final long[] old=set;
						set=new long[old.length*2];
						mask=set.length-1;
						for(long x : old){
							if(x==0){continue;}
							int p=(int)(Tools.hash64shift(x))&mask;
							while(set[p]!=0){p=(p+1)&mask;}
							set[p]=x;
						}
					}
				}
			}
			avll.add(v);
			fll.add(v);
			adds++;
			tr.tick(fll, adds);
			if(distinct>=nextCp){
				report(avll, fll, mle, win, fuse, env, tr, adds, distinct, B, orgs);
				nextCp=Math.max(nextCp*2, distinct+1);
			}
		}
		bf.close();
		report(avll, fll, mle, win, fuse, env, tr, adds, distinct, B, orgs);
		t.stop();
		System.err.println("adds="+adds+" trueDistinct="+distinct+" time="+t);
	}

	static void report(ArithmeticVariableLogLog avll, Fll53 fll, Fll53MLE mle,
			ComplexityFarm.BigNet win, ComplexityFarm.BigNet fuse, double[][] env,
			Fll53Chron.Traj tr, long adds, long distinct, int B, int orgs){
		final double eA=avll.cardinality();
		final double eM=mle.estimate(fll);
		double eB=eM, eC=eM;
		int k=0;
		double ySum=0;
		for(int w0=0; w0+Fll53FarmW.WIN<=orgs; w0+=Fll53FarmW.WIN){
			ySum+=win.predict(Fll53FarmW.features(fll, adds, w0));
			k++;
		}
		if(k>0 && eM>=32.0*B){
			final double cHat=Math.pow(2.0, ySum/k);
			final double w=Math.max(0, Math.min(1, (cHat-0.85)/0.10));
			eB=w*eM+(1-w)*adds*cHat;
			if(fuse!=null){
				final double[] ff=Fll53Chron.fuse(tr, fll, adds, B);
				if(Fll53Chron.inEnvelope(env, ff, 0.05)){
					final double cC=Math.pow(2.0, fuse.predict(ff));
					final double wc=Math.max(0, Math.min(1, (cC-0.85)/0.10));
					eC=wc*eM+(1-wc)*adds*cC;
				}else{eC=eB;}   // envelope guard: out-of-training-range -> blend
			}else{eC=eB;}
		}
		final double d=Math.max(1, distinct);
		System.out.println(String.format(
			"%d\t%d\t%.0f\t%.0f\t%.0f\t%.0f\t%+.2f\t%+.2f\t%+.2f\t%+.2f",
			adds, distinct, eA, eM, eB, eC,
			100*(eA-d)/d, 100*(eM-d)/d, 100*(eB-d)/d, 100*(eC-d)/d));
	}
}
