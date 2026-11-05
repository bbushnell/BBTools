package test.kmer;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;

import dna.AminoAcid;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.HashBuffer;
import kmer.KmerTableSet;
import shared.Shared;
import stream.Streamer;

/**
 * Checks agreement between buffered placement and public table-set count lookup.
 * Literal amino ARA=32 and AAA=0 isolate five-bit core routing; DNA GAA also
 * encodes as32 and supplies the two-bit control. Each fresh JVM uses three fixed
 * 77-cell shards and four-key buffers obtained from an unstarted production loader.
 * No reads are loaded or Java workers started.
 * @author Sucrose
 */
public final class AminoCoreRoutingTest{

	/** Selects the affected buffered case or one of three routing/cell-geometry controls. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=name argument.");
		}
		final String mode=args[0].substring(5);
		if(!mode.equals("amino-core-buffered") && !mode.equals("amino-unmasked-control")
				&& !mode.equals("dna-core-control") && !mode.equals("amino-core-direct-control")){
			throw new IllegalArgumentException("Unknown amino routing case: "+mode);
		}
		final boolean amino=!mode.equals("dna-core-control");
		final boolean core=!mode.equals("amino-unmasked-control");
		final boolean direct=mode.equals("amino-core-direct-control");
		Shared.AMINO_IN=amino;
		AbstractKmerTableSet.MASK_CORE=core;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		AbstractKmerTableSet.DISPLAY_PROGRESS=false;
		check(AminoAcid.acidToNumber['A']==0 && AminoAcid.acidToNumber['R']==1
			&& AminoAcid.baseToNumber['G']==2, mode, "Literal keys require A=0, amino R=1 and DNA G=2.");

		final KmerTableSet tables=new KmerTableSet(new String[] {
			"k=3", "ways=3", "initialsize=17", "threads=1", "buflen=4", "prefilter=0",
			"prealloc=f", "rcomp=f", "minprob=0"}, 12);
		final long queryMask=core ? (amino ? 992L : 12L) : -1L;
		final int queryWay=core && !amino ? 0 : 2;
		check(tables.k==3 && tables.ways==3 && tables.amino==amino && !tables.rcomp()
			&& tables.initialSize==17 && tables.buflen==4 && tables.coreMask==queryMask && tables.middleMask==-1L
			&& tables.cmMask==queryMask && !tables.prefilter && tables.prefilterArray==null
			&& tables.tables()==null, mode, "Fixture requires the requested geometry without automatic allocation.");
		check(tables.kmerToWay(0)==0 && tables.kmerToWay(32)==queryWay,
			mode, "Public query routing must match the literal core-mask geometry.");

		final HashArray1D[] shards=new HashArray1D[3];
		//Match the existing allocator's cell geometry, independently of table-set shard selection.
		for(int i=0; i<shards.length; i++){shards[i]=new HashArray1D(new int[] {17}, core ? 12L : -1L);}
		install(tables, shards);
		check(tables.tables()==shards, mode, "Fixture must install exactly the bounded backend array.");
		final HashBuffer buffer=loaderBuffer(tables, mode);
		if(direct){
			//The two-bit backend cell mask can still retrieve keys placed in the correct query shard.
			check(shards[0].set(0, 1)==1 && shards[2].set(32, 2)==1,
				mode, "Direct control must create the two literal keys in shards0 and2.");
		}else{
			check(buffer.incrementAndReturnNumCreated(32, 1)==0
				&& buffer.incrementAndReturnNumCreated(32, 1)==0
				&& buffer.incrementAndReturnNumCreated(0, 1)==0,
				mode, "Three increments must remain buffered before explicit flush.");
			check(buffer.flush()==2, mode, "Explicit flush must create exactly two distinct keys.");
		}
		verifyStorage(shards, mode);
		check(tables.getCount(0)==1, mode, "Public lookup must retrieve the unaffected AAA count1.");
		final int actual=tables.getCount(32);
		check(actual==2, mode, "Expected count32=2, got="+actual
			+"; writerWay="+buffer.kmerToWay(32)+", queryWay="+tables.kmerToWay(32));
		System.out.println("PASS KMR015 case="+mode);
	}

	/** Installs explicit tiny shards once, solely to avoid the heap-sized production allocator.
	 * Public buffer and getCount methods still perform all tested insertion/query work. */
	private static void install(final KmerTableSet tables, final HashArray1D[] shards){
		try{
			final Field field=KmerTableSet.class.getDeclaredField("tables");
			field.setAccessible(true);
			field.set(tables, shards);
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to install the bounded amino routing fixture.", e);
		}
	}

	/** Uses the actual producer constructor so an amino-aware buffer overload needs no fixture change.
	 * Current construction stores the null stream and allocates only the bounded buffers plus
	 * an unstarted Thread object. Never call start or run: run would dereference the stream. */
	private static HashBuffer loaderBuffer(final KmerTableSet tables, final String mode){
		try{
			final Class<?> type=Class.forName("kmer.KmerTableSet$LoadThread");
			final Constructor<?> constructor=type.getDeclaredConstructor(KmerTableSet.class, Streamer.class);
			constructor.setAccessible(true);
			final Thread loader=(Thread)constructor.newInstance(tables, null);
			check(loader.getState()==Thread.State.NEW && !loader.isAlive(), mode, "Fixture loader must remain unstarted.");
			final Field field=type.getDeclaredField("table");
			field.setAccessible(true);
			return (HashBuffer)field.get(loader);
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to obtain the production loader's bounded buffer.", e);
		}
	}

	/** Requires each literal key exactly once with its count, independently of its physical shard. */
	private static void verifyStorage(final HashArray1D[] shards, final String mode){
		int zeroCopies=0, variantCopies=0;
		long keys=0;
		for(HashArray1D shard : shards){
			check(shard.arrayLength()==77 && shard.victims().size()==0,
				mode, "Fixed shards must not grow or spill into victims.");
			keys+=shard.size();
			if(shard.contains(0)){
				zeroCopies++;
				check(shard.getValue(0)==1, mode, "Stored AAA count must be1.");
			}
			if(shard.contains(32)){
				variantCopies++;
				check(shard.getValue(32)==2, mode, "Stored key32 count must be2.");
			}
		}
		check(keys==2 && zeroCopies==1 && variantCopies==1,
			mode, "Physical storage must contain exactly one copy of each literal key.");
	}

	/** Keeps fixture and count oracles active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR015 case="+mode+": "+message);}
	}
}
