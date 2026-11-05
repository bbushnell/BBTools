package test.kmer;

import java.nio.charset.StandardCharsets;
import java.util.BitSet;

import kmer.HashArray1D;
import kmer.TableLoaderLockFree;
import kmer.TableReader;
import stream.Read;

/**
 * Literal middle-mask checks for matching TableLoaderLockFree/TableReader keys.
 * Loader cases insert two k=5 literals and query both; the isolated reader case
 * manually seeds one masked key without invoking the loader. The short-k control
 * requests middle masking while verifying the supported unmasked short-end path.
 * @author Jean
 */
public final class TableReaderMaskMiddleTest{

	/** Selects the loader/reader path, isolated reader or ordinary unmasked control. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException(
				"Expected case=enabled|reader-enabled|disabled-control|short-mask-control");
		}
		final String mode=args[0].substring(5);
		if(mode.equals("enabled")){enabled();}
		else if(mode.equals("reader-enabled")){readerEnabled();}
		else if(mode.equals("disabled-control")){disabledControl();}
		else if(mode.equals("short-mask-control")){shortMaskControl();}
		else{throw new IllegalArgumentException("Unknown case: "+mode);}
		System.out.println("PASS case="+mode);
	}

	/** Enabling maskMiddle must collapse AAAAA and AACAA onto one encoded key. */
	private static void enabled(){
		final Fixture fixture=new Fixture(true);
		check(fixture.added==1, "KMR111: enabled middle masking inserted "+fixture.added+" keys");
		check(fixture.storedKeys()==1, "Enabled middle masking must store only the masked key");
		check(fixture.value(MASKED_KEY)==1, "Masked key must retain the first literal ID");
		check(fixture.value(UNMASKED_VARIANT_KEY)<1,
			"KMR111: enabled middle masking stored the unmasked variant key");
		check(fixture.query("AAAAA")==1, "The exact first literal must query ID 1");
		check(fixture.query("AACAA")==1, "KMR111: the middle-base variant did not query ID 1");
		fixture.checkFixedShards();
	}

	/** A masked TableReader must query a middle-base variant through the masked key. */
	private static void readerEnabled(){
		final Fixture fixture=new Fixture(true, false);
		check(fixture.value(MASKED_KEY)==1, "Masked key must be manually seeded as ID 1");
		check(fixture.value(UNMASKED_VARIANT_KEY)<1, "The unmasked variant key must be absent");
		check(fixture.query("AAAAA")==1, "The exact first literal must query ID 1");
		check(fixture.query("AACAA")==1, "KMR111: the manually masked variant did not query ID 1");
		fixture.checkFixedShards();
	}

	/** Control: without maskMiddle, the middle-base variant keeps its own key and ID. */
	private static void disabledControl(){
		final Fixture fixture=new Fixture(false);
		check(fixture.added==2, "Disabled middle masking must insert both literal keys");
		check(fixture.storedKeys()==2, "Disabled middle masking must retain two keys");
		check(fixture.value(MASKED_KEY)==1, "AAAAA must retain ID 1 on its unmasked key");
		check(fixture.value(UNMASKED_VARIANT_KEY)==2, "AACAA must retain ID 2 on its unmasked key");
		check(fixture.query("AAAAA")==1, "AAAAA must query ID 1 when masking is disabled");
		check(fixture.query("AACAA")==2, "AACAA must query ID 2 when masking is disabled");
		fixture.checkFixedShards();
	}

	/** Control: short-end mode takes precedence over requested middle masking. */
	private static void shortMaskControl(){
		final Fixture fixture=new Fixture(true, true, SHORT_MINK);
		check(fixture.added==5, "Short-k compatibility must insert five literal keys");
		check(fixture.storedKeys()==5, "Short-k compatibility must retain five keys");
		check(fixture.value(MASKED_KEY)==1, "AAAAA must retain ID 1 on its full key");
		check(fixture.value(UNMASKED_VARIANT_KEY)==2, "AACAA must retain ID 2 on its full key");
		check(fixture.value(AAAA_SHORT_KEY)==1, "AAAA short key must retain ID 1");
		check(fixture.value(AACA_SHORT_KEY)==2, "AACA short key must retain ID 2");
		check(fixture.value(ACAA_SHORT_KEY)==2, "ACAA short key must retain ID 2");
		check(fixture.value(ACCA_SHORT_KEY)<1, "Absent ACCA short key must stay absent");
		check(fixture.query("AAAAA")==1, "AAAAA must query ID 1 in short-k control");
		check(fixture.query("AACAA")==2, "AACAA must query ID 2 in short-k control");
		checkAllMarked(fixture.mark("AACA"), "AACA must mark exactly the short literal");
		check(fixture.mark("ACCA")==null, "Absent ACCA short key must not mark bits");
		fixture.checkFixedShards();
	}

	/** Creates seven tiny shards by loading two literals or seeding one masked key. */
	private static final class Fixture{
		Fixture(final boolean maskMiddle){
			this(maskMiddle, true);
		}

		Fixture(final boolean maskMiddle, final boolean useLoader){
			this(maskMiddle, useLoader, 0);
		}

		Fixture(final boolean maskMiddle, final boolean useLoader, final int mink){
			tables=new HashArray1D[WAYS];
			for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
			if(useLoader){
				final TableLoaderLockFree loader=
					new TableLoaderLockFree(tables, K, mink, 0, 0, 0, false, maskMiddle);
				loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
				added=loader.processData(null, new String[] {"AAAAA", "AACAA"}, true, false, false);
			}else{
				check(tables[(int)(MASKED_KEY%WAYS)].set(MASKED_KEY, 1)==1,
					"Masked key must install as ID 1");
				added=1;
			}
			reader=new TableReader(K, mink, 0, 0, 0, false, maskMiddle);
			reader.countArray=new int[3];
		}

		/** Queries a literal read through the configured TableReader. */
		int query(final String bases){
			return reader.findBestMatch(read(bases), tables);
		}

		/** Marks a read through the short-end-aware masking path. */
		BitSet mark(final String bases){
			return reader.markBits(read(bases), tables);
		}

		/** Looks up one already length-encoded literal key from its fixed shard. */
		int value(final long key){
			return tables[(int)(key%WAYS)].getValue(key);
		}

		/** Counts all keys stored in primary and victim tiers. */
		long storedKeys(){
			long size=0;
			for(HashArray1D table : tables){size+=table.size()+table.victims().size();}
			return size;
		}

		/** Verifies the small shard schedule never resized or spilled while loading. */
		void checkFixedShards(){
			for(HashArray1D table : tables){
				check(table.arrayLength()==77 && table.victims().size()==0,
					"Fixture shards must stay fixed and victim-free");
			}
		}

		final HashArray1D[] tables;
		final TableReader reader;
		final long added;
	}

	/** Creates an unpaired read over the supplied literal bases. */
	private static Read read(final String bases){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, bases, 0);
	}

	/** Keeps literal oracles active under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Requires positions 0..3, and only those positions, to be marked. */
	private static void checkAllMarked(final BitSet bits, final String message){
		check(bits!=null, message);
		check(bits.cardinality()==SHORT_MINK, message);
		for(int i=0; i<SHORT_MINK; i++){check(bits.get(i), message);}
	}

	private static final int K=5;
	private static final int SHORT_MINK=4;
	private static final int WAYS=7;
	private static final long MASKED_KEY=1L<<(2*K);
	private static final long UNMASKED_VARIANT_KEY=MASKED_KEY | (1L<<(2*(K/2)));
	private static final long AAAA_SHORT_KEY=1L<<(2*SHORT_MINK);
	private static final long AACA_SHORT_KEY=AAAA_SHORT_KEY | (1L<<2);
	private static final long ACAA_SHORT_KEY=AAAA_SHORT_KEY | (1L<<(2*(SHORT_MINK/2)));
	private static final long ACCA_SHORT_KEY=AACA_SHORT_KEY | ACAA_SHORT_KEY;
}
