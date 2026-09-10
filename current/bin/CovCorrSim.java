package bin;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.LineParserS1;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Measurement scaffold: synthesizes a correlated-sample cov file from a real
 * independent-sample cov file, mimicking RandomReadsMG's depth model (per-genome,
 * per-sample multiplicative jitter: depth*=exp((U*2-1)*jitter)). Group g of size k
 * emits k near-duplicate columns of source column g. Used to test effective-sample
 * estimators across jitter and group structure without mapping reads
 * (plans/sample_reduction_design_2026-09-09.md, experiment A hardening).
 * @author Amber
 * @date September 9, 2026
 */
public class CovCorrSim {

	public static void main(String[] args){
		Timer t=new Timer();
		CovCorrSim x=new CovCorrSim(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CovCorrSim(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in")){in=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("sizes") || a.equals("groups")){
				String[] parts=b.split(",");
				sizes=new int[parts.length];
				for(int j=0; j<parts.length; j++){sizes[j]=Integer.parseInt(parts[j]);}
			}
			else if(a.equals("jitter")){jitter=Float.parseFloat(b);}
			else if(a.equals("zeroprob")){zeroProb=Float.parseFloat(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		assert(in!=null && out!=null && sizes!=null) : "in=, out=, sizes= required";
	}

	/** Deterministic uniform [0,1) from (seed, tid, outColumn); mirrors RandomReadsMG's
	 * per-file jitter draw being independent per (genome, sample). */
	private double jitterFactor(long tid, int col){
		long h=Tools.hash64shift(seed+0x9E3779B97F4A7C15L*(col+1))^Tools.hash64shift(tid);
		double u=(Tools.hash64shift(h)>>>11)*0x1.0p-53;
		return Math.exp((u*2-1)*jitter);
	}

	/** Absence is a property of (genome, LOGICAL group): duplicates share their zero
	 * pattern; independent groups draw independently. At NEON-like zeroProb (~0.9)
	 * independent samples still share ~zeroProb^2 of rows as zeros — co-absence that
	 * must NOT read as correlation. */
	private boolean absent(long tid, int group){
		if(zeroProb<=0){return false;}
		long h=Tools.hash64shift(seed^0xC2B2AE3D27D4EB4FL)+0x9E3779B97F4A7C15L*(group+1);
		double u=(Tools.hash64shift(h^Tools.hash64shift(tid))>>>11)*0x1.0p-53;
		return u<zeroProb;
	}

	void process(Timer t){
		final int groups=sizes.length;
		int outCols=0;
		for(int s : sizes){outCols+=s;}
		ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, true);
		bsw.start();
		LineParser1 lp=new LineParser1('\t');
		ByteBuilder bb=new ByteBuilder();
		int inDepths=-1;
		long rows=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			bb.clear();
			if(line.length>0 && line[0]=='#'){
				lp.set(line);
				if(lp.termEquals("#Depths", 0)){
					inDepths=lp.parseInt(1);
					assert(inDepths>=groups) : "Source has "+inDepths+" cols; need >= "+groups+" logical groups";
					bb.append("#Depths\t").append(outCols);
				}else if(lp.termEquals("#ShortName", 0)){
					bb.append("#ShortName\tID\tSize");
					for(int c=0; c<outCols; c++){bb.append("\tCov_").append(c);}
					bb.append("\tEdge\tWeight");
				}else{
					bb.append(line);
				}
				bb.nl();
				bsw.print(bb);
				continue;
			}
			lp.set(line);
			String name=lp.parseString(0);
			lpu.set(name);
			//Most contigs are ..._tid_<taxid>; a few are unlabeled (e.g. "contig_104651") and
			//jitter independently under their name hash — correct for an unknown genome.
			final long tid=(lpu.terms()>=2 && lpu.termEquals("tid", lpu.terms()-2)) ?
				lpu.parseLong(lpu.terms()-1) : name.hashCode();
			bb.append(name).tab().append(lp.parseString(1)).tab().append(lp.parseString(2));
			int col=0;
			for(int g=0; g<groups; g++){
				double base=(absent(tid, g) ? 0 : lp.parseFloat(3+g));
				for(int k=0; k<sizes[g]; k++, col++){
					bb.tab().append(base==0 ? 0 : base*jitterFactor(tid, col), 4);
				}
			}
			for(int f=3+inDepths; f<lp.terms(); f++){bb.tab().append(lp.parseString(f));}
			bb.nl();
			bsw.print(bb);
			rows++;
		}
		bf.close();
		bsw.poisonAndWait();
		outstream.println("Wrote "+rows+" rows, "+groups+" groups -> "+outCols+" columns, jitter="+jitter);
		t.stop("Total time:");
	}

	private final LineParserS1 lpu=new LineParserS1('_');
	private String in=null, out=null;
	private int[] sizes=null;
	private float jitter=0.05f;
	private float zeroProb=0f;
	private long seed=1;
	private PrintStream outstream=System.err;
}
