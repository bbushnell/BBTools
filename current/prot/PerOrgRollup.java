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
 * <p>OPTIONAL SPARSE OUTPUT ({@code out2=}, added for Tier B / MANIFEST "TIER_B_DESIGN.md"): a
 * SECOND per-organism file in the 16-field format {@code marker_select.py} (the pre-Java marker
 * precompute) and {@link AntiSetMiner#loadPresence} both parse — {@code tid domain length gc acgt
 * cds mapped glenSum glenSq coding r16 r23 r5 rother trna families}, families as sparse
 * "rank:count;..." (only nonzero ranks), gc/length/etc. as raw SUMMED counts (not the dense path's
 * gc fraction) — verified against both consumers' field indices (domain@1, families@15) before
 * adding this. Written from the SAME per-org accumulation pass as the dense output (one read of
 * the (large) input cache serves both), so out1= behavior is completely unaffected when out2= is
 * omitted — differential-verified byte-identical against the pre-out2 version.
 *
 * <p>Usage: perorgrollup.sh in=percontig_cache.tsv out=perorg_cache.tsv [topn=8000] [out2=perorg_sparse.tsv]
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
			else if(a.equals("out2")){out2=b;}
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

	/** Per-organism accumulator: gene-structure sums + dense family copy-counts (indexed by rank).
	 *  The cds..trna fields exist only to feed the optional sparse (out2=) output — they cost one
	 *  extra int-parse and add per shred regardless, since tracking them unconditionally is simpler
	 *  and cheaper than a branch on whether out2 was requested. */
	static final class Org {
		long length, gc, acgt;
		long cds, mapped, glenSum, glenSq, coding, r16, r23, r5, rother, trna;
		int nshreds;
		String domain;
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

			final String domain=lp.parseString(2);
			if(o.domain==null){o.domain=domain;}
			else{assert(o.domain.equals(domain)) : "Domain changed mid-organism (tid "+tid+"): "
				+o.domain+" -> "+domain+" -- CacheBuilder assigns domain per-tid, this should never vary.";}

			o.length+=lp.parseInt(3);
			o.gc+=lp.parseInt(4);
			o.acgt+=lp.parseInt(5);
			o.cds+=lp.parseInt(6);
			o.mapped+=lp.parseInt(7);
			o.glenSum+=lp.parseLong(8);
			o.glenSq+=lp.parseLong(9);
			o.coding+=lp.parseInt(10);
			o.r16+=lp.parseInt(11);
			o.r23+=lp.parseInt(12);
			o.r5+=lp.parseInt(13);
			o.rother+=lp.parseInt(14);
			o.trna+=lp.parseInt(15);
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

		if(out2!=null){
			FileFormat ffout2=FileFormat.testOutput(out2, FileFormat.TXT, null, true, overwrite, false, false);
			ByteStreamWriter bsw2=new ByteStreamWriter(ffout2);
			bsw2.start();
			int written2=0;
			try{
				ByteBuilder bb=new ByteBuilder(1<<16);
				bb.append("#tid\tdomain\tlength\tgc\tacgt\tcds\tmapped\tglenSum\tglenSq\tcoding"
					+"\tr16\tr23\tr5\trother\ttrna\tfamilies").nl();
				bsw2.print(bb);
				bb.clear();

				for(java.util.Map.Entry<Integer, Org> e : orgs.entrySet()){
					final Org o=e.getValue();
					bb.append(e.getKey().intValue()).append('\t').append(o.domain).append('\t')
						.append(o.length).append('\t').append(o.gc).append('\t').append(o.acgt).append('\t')
						.append(o.cds).append('\t').append(o.mapped).append('\t').append(o.glenSum).append('\t')
						.append(o.glenSq).append('\t').append(o.coding).append('\t').append(o.r16).append('\t')
						.append(o.r23).append('\t').append(o.r5).append('\t').append(o.rother).append('\t')
						.append(o.trna).append('\t');
					boolean first=true;
					for(int i=0; i<topN; i++){
						if(o.fams[i]!=0){
							if(!first){bb.append(';');}
							bb.append(i).append(':').append(o.fams[i]);
							first=false;
						}
					}
					bb.nl();
					bsw2.print(bb);
					bb.clear();
					written2++;
				}
			}finally{
				bsw2.poisonAndWait();
			}
			outstream.println("Wrote "+written2+" organisms (sparse) to "+out2);
		}

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
	private String out2=null;
	private int topN=8000;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
