package jgi;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import assemble.Tadpole;
import kmer.KmerTableSet;

/** Checks ownership budgeting in BBMerge's embedded Tadpole without loading tables.
 * @author Brian Bushnell, Fischl
 */
public final class BBMergeTableMemoryTest {

	/** Creates a tiny input solely for constructor validation, then removes it. */
	public static void main(String[] args) throws Exception{
		assert(args!=null) : "The test launcher must supply an argument array.";
		final Path input=Files.createTempFile("bbmerge-memory-", ".fq");
		try{
			Files.write(input, "@read/1\nACGT\n+\nIIII\n@read/2\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII));
			check(input, 31, false, false);
			check(input, 31, true, false);
			check(input, 31, false, true);
			check(input, 93, false, false);
			check(input, 93, true, false);
			System.out.println("BBMERGE_TABLE_MEMORY_TEST_PASS cases=5");
		}finally{
			Files.delete(input);
		}
	}

	/** Checks the real constructed table owner; reflection avoids a product-only test API. */
	private static void check(final Path input, final int k, final boolean shave,
			final boolean rinse) throws Exception{
		assert(Files.isRegularFile(input)) : "BBMerge constructor requires the generated FASTQ fixture.";
		final BBMerge merge=new BBMerge(new String[] {"in="+input, "out=null", "t=1",
				"net=null", "kfilter=1", "k="+k, "shave="+shave, "rinse="+rinse});
		final Field embedded=BBMerge.class.getDeclaredField("tadpole");
		embedded.setAccessible(true);
		final Tadpole tadpole=(Tadpole)embedded.get(merge);
		final boolean expected=shave || rinse;
		assert(tadpole!=null && tadpole.useOwnership==expected) :
				"BBMerge only allocates ownership during Shaver cleaning; k="+k+
				", shave="+shave+", rinse="+rinse;
		if(k==31){
			final Field bytes=KmerTableSet.class.getDeclaredField("bytesPerKmer");
			bytes.setAccessible(true);
			final int actual=bytes.getInt(tadpole.tables());
			assert(actual==(expected ? 17 : 13)) :
					"Tadpole1 budgets12 key/count bytes,1 existing contig allowance and4 only for ownership; actual="+actual;
		}
	}
}
