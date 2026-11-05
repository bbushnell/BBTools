package test.ukmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import assemble.Contig;
import kmer.AbstractKmerTableSet;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Checks long-key absent coverage in both public overloads with literal k=62 oracles.
 * The only stored key is A^62 with count7; valid absent windows contribute zero but
 * remain in the denominator, while ambiguous windows are skipped.
 * @author Jean */
public final class KmerTableSetUCoverageTest{

	/** Accepts bytes- or contig- followed by mixed, absent, present, ambiguous, no-valid or short. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=<bytes|contig>-<scenario> argument.");
		}
		final String mode=args[0].substring(5);
		final boolean contig=mode.startsWith("contig-");
		if(!contig && !mode.startsWith("bytes-")){throw new IllegalArgumentException("Unknown coverage-long case: "+mode);}
		final String scenario=mode.substring(contig ? 7 : 6);
		final String sequence;
		final float average;
		final int minimum, maximum;
		switch(scenario){
			case "mixed": sequence=repeat('A', K)+"C"; average=3.5f; minimum=0; maximum=7; break;
			case "absent": sequence=repeat('C', K); average=0; minimum=0; maximum=0; break;
			case "present": sequence=repeat('A', K+1); average=7; minimum=7; maximum=7; break;
			case "ambiguous": sequence=repeat('A', K)+"N"+repeat('C', K); average=3.5f; minimum=0; maximum=7; break;
			case "no-valid": sequence=repeat('N', K); average=0; minimum=0; maximum=0; break;
			case "short": sequence=repeat('A', K-1); average=0; minimum=0; maximum=0; break;
			default: throw new IllegalArgumentException("Unknown coverage-long scenario: "+scenario);
		}

		final Field orientation=Kmer.class.getDeclaredField("rcomp");
		orientation.setAccessible(true);
		final boolean oldPacked=Kmer.PACKED, oldRcomp=orientation.getBoolean(null);
		final boolean oldMaskCore=AbstractKmerTableSet.MASK_CORE;
		final boolean oldMaskMiddle=AbstractKmerTableSet.MASK_MIDDLE;
		try{
			Kmer.PACKED=false;
			AbstractKmerTableSet.MASK_CORE=false;
			AbstractKmerTableSet.MASK_MIDDLE=false;
			final FixtureSet set=new FixtureSet("in=/dev/null", "k=62", "ways=1", "threads=1",
				"prealloc=f", "prefilter=0", "onepass=f", "rcomp=f", "qin=33", "interleaved=f");
			final HashArrayU1D table=new HashArrayU1D(new int[] {101}, set.k, set.kbig);
			set.installForTest(table);
			final Kmer present=kmer(set, repeat('A', K), mode);
			coverage(table.set(present, COUNT)==1, mode, "Stored A^62 key must be new");
			verifyTables(set, table, mode);

			if(contig){
				checkContig(set, sequence, average, minimum, maximum, scenario.equals("short"), mode);
				if(scenario.equals("short")){checkContig(set, "", 0, 0, 0, true, mode);}
			}else{
				checkBytes(set, sequence, average, mode);
				if(scenario.equals("short")){checkBytes(set, "", 0, mode);}
			}
			verifyTables(set, table, mode);
		}finally{
			AbstractKmerTableSet.MASK_MIDDLE=oldMaskMiddle;
			AbstractKmerTableSet.MASK_CORE=oldMaskCore;
			orientation.setBoolean(null, oldRcomp);
			Kmer.PACKED=oldPacked;
		}
		System.out.println("PASS coverage-long case="+mode+" average="+average);
	}

	/** Checks the byte-array overload and unchanged input bytes. */
	private static void checkBytes(final KmerTableSetU set, final String sequence,
			final float expected, final String mode){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] before=bases.clone();
		final float actual=set.calcCoverage(bases, bases.length, new Kmer(set.kbig()));
		coverage(actual==expected, mode, "Expected average="+expected+", actual="+actual);
		coverage(Arrays.equals(bases, before), mode, "Coverage must not change the query bases");
	}

	/** Checks returned and stored statistics, including the existing short-contig early return. */
	private static void checkContig(final KmerTableSetU set, final String sequence,
			final float expected, final int minimum, final int maximum, final boolean shortInput,
			final String mode){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] before=bases.clone();
		final Contig contig=new Contig(bases);
		contig.coverage=9;
		contig.minCov=7;
		contig.maxCov=11;
		final float actual=set.calcCoverage(contig, new Kmer(set.kbig()));
		coverage(actual==expected, mode, "Expected average="+expected+", actual="+actual);
		if(shortInput){
			coverage(contig.coverage==9 && contig.minCov==7 && contig.maxCov==11, mode,
				"Short contig must return zero without overwriting existing statistics");
		}else{
			coverage(contig.coverage==expected && contig.minCov==minimum && contig.maxCov==maximum, mode,
				"Expected coverage/min/max="+expected+"/"+minimum+"/"+maximum+", actual="+
				contig.coverage+"/"+contig.minCov+"/"+contig.maxCov);
		}
		coverage(contig.bases==bases && Arrays.equals(bases, before), mode,
			"Coverage must preserve the contig's base array and its contents");
	}

	/** Guards actual tiny storage and keeps the raw lookup absence sentinel unchanged. */
	private static void verifyTables(final KmerTableSetU set, final HashArrayU1D table,
			final String mode){
		coverage(set.kbig()==K && set.tables().length==1 && set.tables()[0]==table
			&& table.arrayLength()==161 && table.victims().size()==0 && table.size()==1, mode,
			"Coverage-long fixture must retain one fixed161-cell shard and one key without victims");
		coverage(set.getCount(kmer(set, repeat('A', K), mode))==COUNT
			&& set.getCount(kmer(set, repeat('C', K), mode))==-1
			&& set.getCount(kmer(set, repeat('A', K-1)+"C", mode))==-1, mode,
			"A^62 must keep count7 and absent long keys must still return -1");
	}

	/** Installs a bounded table shard without invoking ScheduleMaker allocation. */
	private static final class FixtureSet extends KmerTableSetU{
		FixtureSet(String... args){super(args, 0);}

		void installForTest(final AbstractKmerTableU table) throws Exception{
			TABLES.set(this, new AbstractKmerTableU[] {table});
		}

		private static final Field TABLES;
		static{
			try{
				TABLES=KmerTableSetU.class.getDeclaredField("tables");
				TABLES.setAccessible(true);
			}catch(Exception e){throw new ExceptionInInitializerError(e);}
		}
	}

	/** Builds one complete long k-mer through the production rolling API. */
	private static Kmer kmer(final KmerTableSetU set, final String sequence, final String mode){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final Kmer kmer=new Kmer(set.kbig());
		for(int i=0; i<bases.length; i++){kmer.addRight(bases[i]);}
		coverage(kmer.len()==set.kbig(), mode, "Fixture sequence must build one full k-mer: "+sequence.length());
		return kmer;
	}

	/** Creates a literal sequence without relying on Java 11 String.repeat. */
	private static String repeat(final char c, final int count){
		final char[] chars=new char[count];
		Arrays.fill(chars, c);
		return new String(chars);
	}

	/** Keeps literal coverage expectations active with production assertions disabled. */
	private static void coverage(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("UKM113 case="+mode+": "+message);}
	}

	private static final int K=62;
	private static final int COUNT=7;
}
