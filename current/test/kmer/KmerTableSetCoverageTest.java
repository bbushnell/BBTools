package test.kmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import assemble.Contig;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.KmerTableSet;
import shared.Shared;

/** Checks absent-key coverage in both public overloads with literal k=3 oracles.
 * AAA has count4 and AAT count8; AAC and ACG are absent. Valid absent windows
 * contribute zero but remain in the denominator; ambiguous windows are skipped.
 * Uses the existing explicit-table fixture pattern without input loading or
 * automatic allocation. Run one case per fresh JVM under either assertion mode.
 * @author Sucrose
 */
public final class KmerTableSetCoverageTest{

	/** Accepts bytes- or contig- followed by mixed, absent, present, ambiguous, zero-sum, no-valid or short. */
	public static void main(String[] args) throws ReflectiveOperationException{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=<bytes|contig>-<scenario> argument.");
		}
		final String mode=args[0].substring(5);
		final boolean contig=mode.startsWith("contig-");
		if(!contig && !mode.startsWith("bytes-")){throw new IllegalArgumentException("Unknown coverage case: "+mode);}
		final String scenario=mode.substring(contig ? 7 : 6);
		final String sequence;
		final float average;
		final int minimum, maximum;
		switch(scenario){
			case "mixed": sequence="AAAC"; average=2; minimum=0; maximum=4; break;
			case "absent": sequence="ACG"; average=0; minimum=0; maximum=0; break;
			case "present": sequence="AAAT"; average=6; minimum=4; maximum=8; break;
			case "ambiguous": sequence="AAANAAC"; average=2; minimum=0; maximum=4; break;
			case "zero-sum": sequence="AAACGCA"; average=0.8f; minimum=0; maximum=4; break;
			case "no-valid": sequence="ANAN"; average=0; minimum=0; maximum=0; break;
			case "short": sequence="AA"; average=0; minimum=0; maximum=0; break;
			default: throw new IllegalArgumentException("Unknown coverage scenario: "+scenario);
		}

		final boolean oldAmino=Shared.AMINO_IN, oldCore=AbstractKmerTableSet.MASK_CORE;
		final boolean oldMiddle=AbstractKmerTableSet.MASK_MIDDLE, oldProgress=AbstractKmerTableSet.DISPLAY_PROGRESS;
		try{
			Shared.AMINO_IN=false;
			AbstractKmerTableSet.MASK_CORE=false;
			AbstractKmerTableSet.MASK_MIDDLE=false;
			AbstractKmerTableSet.DISPLAY_PROGRESS=false;
			final KmerTableSet set=new KmerTableSet(new String[] {
				"k=3", "ways=1", "initialsize=17", "prefilter=0", "onepass=f",
				"prealloc=f", "rcomp=t", "minprob=0"}, 12);
			fixture(set.k==3 && set.ways==1 && !set.amino && set.rcomp() && !set.prefilter
				&& !set.onePass && set.tables()==null && set.prefilterArray==null
				&& set.in1.isEmpty() && set.in2.isEmpty() && set.extra.isEmpty(), mode,
				"Fixture requires unallocated canonical DNA tables without inputs or prefilter");
			final HashArray1D table=new HashArray1D(new int[] {17}, -1L);
			//AAA: forward0, reverse63; AAT: forward3, reverse15. TableSet keys have no length tag.
			fixture(table.set(63, 4)==1 && table.set(15, 8)==1, mode, "Both literal count keys must be newly stored");
			final AbstractKmerTable[] shards={table};
			install(set, shards);
			verifyTables(set, shards, table, mode);

			if(contig){
				checkContig(set, sequence, average, minimum, maximum, scenario.equals("short"), mode);
				if(scenario.equals("short")){checkContig(set, "", 0, 0, 0, true, mode);}
			}else{
				//The present control also verifies that blength bounds the scanned prefix.
				final String backing=scenario.equals("present") ? "AAATACG" : sequence;
				checkBytes(set, backing, sequence.length(), average, mode);
				if(scenario.equals("short")){checkBytes(set, "", 0, 0, mode);}
			}
			verifyTables(set, shards, table, mode);
		}finally{
			Shared.AMINO_IN=oldAmino;
			AbstractKmerTableSet.MASK_CORE=oldCore;
			AbstractKmerTableSet.MASK_MIDDLE=oldMiddle;
			AbstractKmerTableSet.DISPLAY_PROGRESS=oldProgress;
		}
		System.out.println("PASS coverage case="+mode+" average="+average);
	}

	/** Installs one explicitly sized shard, as in the existing TableSet worker fixture. */
	private static void install(final KmerTableSet set, final AbstractKmerTable[] shards)
			throws ReflectiveOperationException{
		final Field field=KmerTableSet.class.getDeclaredField("tables");
		field.setAccessible(true);
		field.set(set, shards);
	}

	/** Checks the byte-array overload, bounded prefix and unchanged input bytes. */
	private static void checkBytes(final KmerTableSet set, final String sequence,
			final int length, final float expected, final String mode){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] before=bases.clone();
		final float actual=set.calcCoverage(bases, length);
		coverage(actual==expected, mode, "Expected average="+expected+", actual="+actual);
		fixture(Arrays.equals(bases, before), mode, "Coverage must not change the query bases");
	}

	/** Checks returned and stored statistics, including the existing short-contig early return. */
	private static void checkContig(final KmerTableSet set, final String sequence,
			final float expected, final int minimum, final int maximum, final boolean shortInput,
			final String mode){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] before=bases.clone();
		final Contig contig=new Contig(bases);
		contig.coverage=9;
		contig.minCov=7;
		contig.maxCov=11;
		final float actual=set.calcCoverage(contig);
		coverage(actual==expected, mode, "Expected average="+expected+", actual="+actual);
		if(shortInput){
			coverage(contig.coverage==9 && contig.minCov==7 && contig.maxCov==11, mode,
				"Short contig must return zero without overwriting existing statistics");
		}else{
			coverage(contig.coverage==expected && contig.minCov==minimum && contig.maxCov==maximum, mode,
				"Expected coverage/min/max="+expected+"/"+minimum+"/"+maximum+", actual="+
				contig.coverage+"/"+contig.minCov+"/"+contig.maxCov);
		}
		fixture(contig.bases==bases && Arrays.equals(bases, before), mode,
			"Coverage must preserve the contig's base array and its contents");
	}

	/** Guards actual tiny storage and keeps the raw lookup absence sentinel unchanged. */
	private static void verifyTables(final KmerTableSet set, final AbstractKmerTable[] shards,
			final HashArray1D table, final String mode){
		fixture(set.tables()==shards && shards.length==1 && shards[0]==table
			&& table.arrayLength()==77 && table.victims().size()==0 && table.size()==2, mode,
			"Coverage fixture must retain one fixed77-cell shard and exactly two keys without victims");
		fixture(table.getValue(63)==4 && table.getValue(15)==8
			&& set.getCount(0, 63)==4 && set.getCount(3, 15)==8, mode,
			"Literal AAA/AAT keys must keep counts4/8 under the real canonical lookup");
		fixture(set.getCount(1, 47)==-1 && set.getCount(6, 27)==-1
			&& set.getCount(25, 38)==-1 && set.getCount(36, 57)==-1, mode,
			"Raw AAC/ACG/CGC/GCA lookup must still return -1; only coverage clamps absent counts");
		fixture(!set.errorState && set.readsIn==0 && set.basesIn==0 && set.kmersIn==0
			&& set.prefilterArray==null, mode, "Direct coverage must not load data or create a prefilter");
	}

	/** Keeps literal coverage expectations active with production assertions disabled. */
	private static void coverage(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR022 case="+mode+": "+message);}
	}

	/** Distinguishes fixture setup/storage failures from the intended missing-key discrepancy. */
	private static void fixture(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("Fixture qualification case="+mode+": "+message);}
	}
}
