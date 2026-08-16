package prot;

import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Aggregates the 17-field per-shred cache ({@link CacheBuilder}'s output) into a per-organism
 * cache by summing gene-structure stats and merging sparse family copy-counts across every shred
 * sharing a taxonomic ID.
 *
 * <p>OUTPUT SCHEMA IS UNCHANGED from the pre-rewrite version — this is an input-format migration
 * only, not a new contract: {@code #tid length gc nshreds fam_0..fam_(topN-1)}, dense. Only the
 * INPUT side changed (17-field sparse-family cache, not the dead 4+dense Python format), per
 * CACHEBUILDER_DESIGN.md §7. The dead {@code tid=-1} accommodation (tidmap fallback) is DROPPED
 * entirely — the corpus has carried zero {@code tid=-1} shreds since the 08 source fix, and a
 * non-positive tid here is now a real anomaly, not an expected case.
 *
 * <p>gc column = Sum(shred gc-count) / Sum(shred acgt-count), a %.4f fraction — the correct
 * weighted average given the new cache's gc/acgt are raw INT COUNTS, not a pre-computed per-shred
 * fraction (which is what the old dense-format input carried and length-weighted-averaged).
 *
 * <p>Usage: perorgrollup.sh in=percontig_cache.tsv out=perorg_cache.tsv [topn=8000]
 *
 * @author Brian Bushnell, Eru
 */
public class PerOrgRollup {

	public static void main(String[] args){
		Timer t=new Timer();
		PerOrgRollup x=new PerOrgRollup(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public PerOrgRollup(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("topn") || a.equals("top")){topN=Integer.parseInt(b);}
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		{
			Parser.processQuality();
			overwrite=parser.overwrite;
			in1=parser.in1;
			out1=parser.out1;
		}
		assert(in1!=null) : "No input file specified.";
		assert(out1!=null) : "No output file specified.";
	}

	/** Per-organism accumulator: gene-structure sums + dense family copy-counts (indexed by rank). */
	static final class Org {
		long length, gc, acgt;
		int nshreds;
		final long[] fams;
		Org(int topN){fams=new long[topN];}
	}

	void process(Timer t){
		final TreeMap<Integer, Org> orgs=new TreeMap<Integer, Org>();
		final ByteFile bf=ByteFile.makeByteFile(in1, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final IntList rankBuf=new IntList(64), countBuf=new IntList(64);

		long rowsRead=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			assert(lp.terms()>=17) : "Malformed cache row: "+lp.terms()+" fields (need >=17): "+new String(line);

			final int tid=lp.parseInt(1);
			assert(tid>0) : "Non-positive tid in cache row (corpus should carry only valid tids "
				+"after the source fix): "+new String(line);

			Org o=orgs.get(tid);
			if(o==null){o=new Org(topN); orgs.put(tid, o);}

			o.length+=lp.parseInt(3);
			o.gc+=lp.parseInt(4);
			o.acgt+=lp.parseInt(5);
			o.nshreds++;

			final int flen=lp.length(16);
			if(flen>0){
				rankBuf.clear(); countBuf.clear();
				parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
				for(int i=0; i<rankBuf.size(); i++){
					final int rank=rankBuf.array[i];
					if(rank<topN){o.fams[rank]+=countBuf.array[i];}
				}
			}

			rowsRead++;
			if(rowsRead%500000==0){outstream.println("  "+rowsRead+" rows read, "+orgs.size()+" orgs");}
		}
		bf.close();

		outstream.println("Rows read: "+rowsRead);
		outstream.println("Orgs: "+orgs.size());

		FileFormat ffout=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		int written=0;
		try{
			ByteBuilder bb=new ByteBuilder(1<<16);
			bb.append("#tid\tlength\tgc\tnshreds");
			for(int i=0; i<topN; i++){bb.append('\t').append("fam_").append(i);}
			bb.nl();
			bsw.print(bb);
			bb.clear();

			for(java.util.Map.Entry<Integer, Org> e : orgs.entrySet()){
				final Org o=e.getValue();
				final double gcFrac=o.acgt>0 ? (double)o.gc/o.acgt : 0;

				bb.append(e.getKey().intValue()).append('\t').append(o.length).append('\t');
				bb.append(gcFrac, 4).append('\t').append(o.nshreds);
				for(int i=0; i<topN; i++){
					bb.append('\t').append(o.fams[i]);
				}
				bb.nl();
				bsw.print(bb);
				bb.clear();
				written++;
			}
		}finally{
			bsw.poisonAndWait();
		}

		outstream.println("Wrote "+written+" organisms to "+out1);
		t.stop();
		outstream.println("Time: \t"+t);
	}

	/** Parses "rank:count;rank:count;..." within [a,b) of line into the (cleared) output
	 *  lists, in order. Zero allocation beyond the two IntList's own backing-array growth.
	 *  Mirrors MagQCVectorMaker.parseFamCounts exactly (same cache field, same encoding). */
	private static void parseFamCounts(byte[] line, int a, int b, IntList outRank, IntList outCount){
		int start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==';'){
				if(i>start){
					int colon=-1;
					for(int j=start; j<i; j++){if(line[j]==':'){colon=j; break;}}
					assert(colon>start) : "Malformed famcounts field (missing ':' in a rank:count pair).";
					outRank.add(parse.Parse.parseInt(line, start, colon));
					outCount.add(parse.Parse.parseInt(line, colon+1, i));
				}
				start=i+1;
			}
		}
	}

	private String in1=null;
	private String out1=null;
	private int topN=8000;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
