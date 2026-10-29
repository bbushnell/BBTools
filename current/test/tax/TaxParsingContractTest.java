package test.tax;

import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.Lineage;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite numeric-header and lineage parsing/reset checks with literal expected results.
 * Each invocation owns its JVM's shared taxonomy state and uses only a tiny local tree.
 * @author Brian Bushnell, Shinobu */
public final class TaxParsingContractTest{

	/** Runs baseline or fixed expectations alongside unchanged parsing controls.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If the temporary tree or state inspection fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		checkNumericBoundaries(baseline);
		check(TaxTree.parseTaxID("tid|4|First")==4, "Literal tid header control failed");
		check(TaxTree.parseTaxID("prefix_ncbi_5.fna")==5, "Filename header control failed");
		check(TaxTree.parseTaxID("plastid_name_tid_4")==4, "False substring must not hide a later valid ID");
		check(TaxTree.parseTaxID("tid|2147483648")<0, "An out-of-int ID must not be accepted");
		// This decimal is 2^64+2; an unchecked long accumulator wraps to taxID2.
		final int overflow=TaxTree.parseTaxID("tid|18446744073709551618");
		check(overflow==(baseline ? 2 : -1), "Overflowed header ID: "+overflow);
		final int later=TaxTree.parseTaxID("tid|18446744073709551618_tid|4");
		check(later==(baseline ? 2 : 4), "Invalid numeric token must not hide the later valid token: "+later);

		check(TaxTree.getTree()==null, "Run this fixture in its own JVM; it owns the shared tree");
		final Path file=Files.createTempFile("tax-parsing-contract-", ".tsv");
		boolean shortFailed=false;
		try{
			TaxTreeText.write(TaxTreeContractTest.tree(TaxTree.GENUS_E), file.toString(), true);
			TaxTree.ALLOW_SHARED_TREE=true;
			TaxTree.LOAD_FAKE_TREE=false;
			final TaxTree tree=TaxTree.loadTaxTree(file.toString(), null, true, false);
			check(tree.getNode(4).name.equals("First"), "Literal fixture leaf failed to load");
			final Lineage full=new Lineage("species__First");
			check(full.nodes[TaxTree.SPECIES]==tree.getNode(4), "Existing full-rank spelling must still resolve");
			try{
				final Lineage abbreviated=new Lineage("s__First");
				check(abbreviated.nodes[TaxTree.SPECIES]==tree.getNode(4), "Short rank must resolve the same species");
			}catch(NullPointerException e){
				if(!baseline){throw e;}
				shortFailed=true;
			}
			check(shortFailed==baseline, "Short-rank lookup should fail only in the baseline");
			final Lineage reset=new Lineage(4);
			check(level(reset)==TaxTree.SPECIES, "Resolved lineage must carry its node's normal rank");
			reset.clear();
			check(reset.line==null, "Clear must discard original lineage text");
			for(Object node : reset.nodes){check(node==null, "Clear must discard all node references");}
			check(level(reset)==(baseline ? TaxTree.SPECIES : -1), "Clear retained a stale rank");
			reset.set(4);
			reset.set(0);//The small fixture has no node in its zero slot.
			check(level(reset)==(baseline ? TaxTree.SPECIES : -1), "An unresolved ID retained a stale rank");
		}finally{Files.deleteIfExists(file);}
		System.out.println("PASS tax_parsing baseline="+baseline+" overflow="+overflow+
				" later="+later+" short_rank_failed="+shortFailed);
	}

	/** Prevents construction of this standalone fixture. */
	private TaxParsingContractTest(){}

	/** Checks literal numeric boundaries and overflow for every supported marker.
	 * @param baseline Expect unchecked long wraparound only in the original parser */
	private static void checkNumericBoundaries(boolean baseline){
		for(String prefix : new String[]{"tid|", "tid_", "ncbi|", "ncbi_"}){
			check(TaxTree.parseTaxID(prefix+"000000000000000000000000000004")==4,
					"Leading-zero control failed for "+prefix);
			check(TaxTree.parseTaxID(prefix+"2147483646")==2147483646,
					"Largest accepted ID changed for "+prefix);
			check(TaxTree.parseTaxID(prefix+"00000000002147483646")==2147483646,
					"Leading zeroes must preserve the largest accepted ID for "+prefix);
			check(TaxTree.parseTaxID(prefix+"2147483647")==-1,
					"Reserved int maximum must remain rejected for "+prefix);
			check(TaxTree.parseTaxID(prefix+"00000000002147483647")==-1,
					"Leading zeroes must not admit the reserved maximum for "+prefix);
			check(TaxTree.parseTaxID(prefix+"0")==-1,
					"Zero must remain rejected for "+prefix);
			check(TaxTree.parseTaxID(prefix+"000000000000000000000000")==-1,
					"An all-zero token must remain rejected for "+prefix);
			check(TaxTree.parseTaxID(prefix+"000000000000000000000000_"+prefix+"4")==4,
					"An all-zero token must not hide the later valid marker for "+prefix);
			check(TaxTree.parseTaxID(prefix+"4294967296")==-1,
					"A value with bit32 set must remain rejected for "+prefix);
			check(TaxTree.parseTaxID(prefix+"18446744073709551617")==(baseline ? 1 : -1),
					"Wrapped 2^64+1 must be rejected after repair for "+prefix);
			check(TaxTree.parseTaxID(prefix+"18446744073709551617_"+prefix+"4")==(baseline ? 1 : 4),
					"An invalid long token must not hide the later marker for "+prefix);
		}
	}

	/** Reads the existing package-private state without widening the production API.
	 * @param lineage View whose retained rank is being checked
	 * @return Stored normal rank
	 * @throws ReflectiveOperationException If the documented state field is unavailable */
	private static int level(Lineage lineage) throws ReflectiveOperationException{
		final Field field=Lineage.class.getDeclaredField("level");
		field.setAccessible(true);
		return field.getInt(lineage);
	}

	/** Keeps independent checks active under -da.
	 * @param condition Required result
	 * @param message Failure explanation */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
