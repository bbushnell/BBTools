package prok;

import java.util.Locale;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.LongIntMap;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Exact sense-oriented member coverage at a declared distinct17mer minimum.
 * Preserves every input copy; no window support, reverse complements or dedupe.
 * @author Raiden
 */
public final class Euk18sSeedMinimumCoverage {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,Euk18sSeedMinimumCoverage.class,false);
		String in=null,seeds=null,clade=null,out=null;int minimum=-1;long expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}else if(key.equals("seeds")){seeds=value;}
			else if(key.equals("clade")){clade=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("minhits")){minimum=Parse.parseIntKMG(value);}
			else if(key.equals("expected")){expected=Parse.parseKMG(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(in!=null && seeds!=null && clade!=null && out!=null && minimum>0 && expected>0,
				"Declare exact-member FASTA, seed FASTA, clade, fresh output, positive minhits and expected count");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		new Euk18sSeedMinimumCoverage(ProkObject.loadLongKmers(seeds,17)).run(in,clade,out,minimum,expected);
		Shared.closeStream(pp.outstream);
	}

	public Euk18sSeedMinimumCoverage(LongHashSet keys){
		require(keys!=null && keys.size()>0,"Coverage must use a nonempty actual caller-visible key set");seeds=keys;
	}

	public void run(String fasta,String clade,String prefix,int minimum,long expected){
		require(minimum>0 && expected>0 && clade!=null && !clade.isEmpty() && clade.indexOf('\t')<0,
				"Coverage denominator and minimum must be explicit; clade is one TSV field");
		RrnaCmWindowParity.fresh(prefix,".summary.tsv",".below.tsv");
		final long[][] counts=new long[LENGTHS.length][2];
		final ByteStreamWriter below=RrnaSourceTargetPrep.writer(prefix+".below.tsv");
		final ByteBuilder b=new ByteBuilder();final long[] rows={0,0};
		try{
			below.println("clade\tinputOrdinal\tassembly\tmember\tlength\tlengthClass\tdistinctKeys\tminimum");
			RrnaResourceIO.stream(fasta,r->{
				final int hits=count(r.bases),group=lengthClass(r.length()),covered=hits>=minimum?1:0;
				counts[0][covered]++;counts[group][covered]++;rows[0]++;
				if(covered==0){
					final String assembly=assemblyKey(r.id);
					below.print(b.clear().append(clade).tab().append(rows[0]).tab().append(assembly).tab()
							.append(RrnaCallerModelAssignment.first(r.id)).tab().append(r.length()).tab()
							.append(LENGTHS[group]).tab().append(hits).tab().append(minimum).nl());rows[1]++;
				}
			});
		}finally{RrnaCmWindowParity.finish(below);}
		require(rows[0]==expected && counts[0][0]+counts[0][1]==expected && rows[1]==counts[0][0],
				"Every raw copy must occur once in covered/below-min accounting; observed="+rows[0]+", expected="+expected);
		long lower=0,upper=0;for(int i=1;i<counts.length;i++){lower+=counts[i][0];upper+=counts[i][1];}
		require(lower==counts[0][0] && upper==counts[0][1],"Length bins must partition both below-min and covered members");
		final ByteStreamWriter summary=RrnaSourceTargetPrep.writer(prefix+".summary.tsv");
		try{
			summary.println("clade\tlengthClass\tminimum\tloadedKeys\tmembers\tcovered\tbelowMinimum");
			for(int i=0;i<counts.length;i++){
				summary.print(b.clear().append(clade).tab().append(LENGTHS[i]).tab().append(minimum).tab()
						.append(seeds.size()).tab().append(counts[i][0]+counts[i][1]).tab()
						.append(counts[i][1]).tab().append(counts[i][0]).nl());
			}
		}finally{RrnaCmWindowParity.finish(summary);}
		System.err.println("EUK18S_SEED_MINIMUM_PASS clade="+clade+" members="+expected+" minimum="+minimum+
				" loadedKeys="+seeds.size()+" covered="+counts[0][1]+" below="+counts[0][0]);
	}

	public int count(byte[] bases){
		return count(bases, null);
	}
	/** Optional member-support audit: one increment per distinct matched key in
	 * this member, irrespective of repeated occurrences. The caller owns the map. */
	public int count(byte[] bases, LongIntMap memberSupport){
		require(bases!=null && bases.length>0,"Exact-gene coverage cannot silently count a missing sequence");
		matched.clear();long word=0;int valid=0;
		for(byte base:bases){
			final int code=AminoAcid.baseToNumber[base&255];
			if(code<0){word=0;valid=0;continue;}
			word=((word<<2)|code)&MASK;
			if(++valid>=17 && seeds.contains(word) && matched.add(word) && memberSupport!=null){memberSupport.increment(word);}
		}
		assert(matched.size()<=seeds.size() && matched.size()<=Math.max(0,bases.length-16)):
				"Only loaded keys at valid in-gene17mer positions may contribute distinct support";
		return matched.size();
	}

	public static int lengthClass(int length){
		require(length>0,"Length-stratified member coverage requires a nonempty exact gene");
		return length<916?1:length<1600?2:length<=2200?3:length<=2500?4:5;
	}
	/** Preserve the existing space-delimited assembly/genome convention without importing unrelated analyses. */
	private static String assemblyKey(final String header){
		assert(header!=null) : "Below-minimum output rows require the original FASTA header metadata";
		int at=header.indexOf(" assembly="), width=10;
		if(at<0){at=header.indexOf(" genome=");width=8;}
		require(at>=0,"Exact-gene FASTA requires explicit assembly= or genome= metadata: "+header);
		final int start=at+width;int stop=start;
		while(stop<header.length() && !Character.isWhitespace(header.charAt(stop))){stop++;}
		require(stop>start,"Assembly metadata must be nonempty: "+header);
		return header.substring(start,stop);
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String[] LENGTHS={"ALL","below916","916_1599","1600_2200","2201_2500","over2500"};
	static final long MASK=(1L<<34)-1;
	final LongHashSet seeds,matched=new LongHashSet(128);
}
