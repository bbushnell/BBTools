package test.tax;

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;

import tax.TaxFilter;
import tax.TaxTree;

/** Finite baseline for TaxFilter copy, level-lowering and substring matching contracts.
 * @author Brian Bushnell, Sucrose */
public final class TaxFilterContractTest{

	/** Runs baseline or fixed expectations.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If the tiny TaxTree fixture cannot be built */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final boolean assertions=assertionsEnabled();
		final TaxTree tree=TaxTreeContractTest.tree(TaxTree.GENUS_E);

		checkDeepCopy(tree, baseline);
		checkLevelLowering(tree, baseline, assertions);
		checkMixedCaseContains(tree, baseline);
		checkTextFileNormalReturnLeaks(tree, baseline);
		if(!baseline){checkAddNumbersFileParseFailureCloses(tree);}

		System.out.println("TAXFILTER_CONTRACT_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed")+" assertions="+assertions);
	}

	/** Prevents construction of this standalone fixture. */
	private TaxFilterContractTest(){}

	/** Verifies that deepCopy returns an independent copy of the selected taxID set.
	 * @param tree Tiny taxonomy tree
	 * @param baseline True to expect the non-Cloneable baseline failure */
	private static void checkDeepCopy(final TaxTree tree, final boolean baseline){
		final TaxFilter filter=new TaxFilter(tree, true);
		filter.addNumber(4, false);
		final TaxFilter copy=filter.deepCopy();
		if(baseline){
			check(copy==null, "Baseline TaxFilter.deepCopy must return null");
		}else{
			check(copy!=null, "Fixed TaxFilter.deepCopy must return a copy");
			check(copy!=filter, "Deep copy must return a distinct filter");
			check(copy.tree()==tree, "Deep copy must share the intended tree");
			checkTaxIDs(copy, 4);
			filter.addNumber(5, false);
			checkTaxIDs(filter, 4, 5);
			checkTaxIDs(copy, 4);
			copy.addNumber(2, false);
			checkTaxIDs(filter, 4, 5);
			checkTaxIDs(copy, 2, 4);
		}
	}

	/** Verifies that lowering the level is allowed for empty nonnull and null tax sets.
	 * @param tree Tiny taxonomy tree
	 * @param baseline True to expect the assertion-only baseline failures
	 * @param assertions Whether this JVM evaluates assertions */
	private static void checkLevelLowering(final TaxTree tree, final boolean baseline,
			final boolean assertions){
		final Throwable emptyFailure=lowerEmptySet(tree, false);
		final Throwable nullFailure=lowerEmptySet(tree, true);
		if(baseline && assertions){
			check(emptyFailure instanceof AssertionError,
					"Baseline empty-set lowering should trip an assertion: "+emptyFailure);
			check(nullFailure instanceof NullPointerException,
					"Baseline null-set lowering should dereference null: "+nullFailure);
		}else{
			check(emptyFailure==null, "Empty-set level lowering should succeed: "+emptyFailure);
			check(nullFailure==null, "Null-set level lowering should succeed: "+nullFailure);
		}
	}

	/** Attempts one level-lowering operation.
	 * @param tree Tiny taxonomy tree
	 * @param nullSet True to clear the selected taxID set before lowering
	 * @return Thrown failure, or null on success */
	private static Throwable lowerEmptySet(final TaxTree tree, final boolean nullSet){
		final TaxFilter filter=new TaxFilter(tree, true);
		filter.setLevel(TaxTree.GENUS, false);
		if(nullSet){filter.clearSet();}
		try{
			filter.setLevel(TaxTree.SPECIES, false);
			check(filter.taxLevel()==TaxTree.SPECIES,
					"Level lowering must set species level: "+filter.taxLevel());
			return null;
		}catch(Throwable t){return t;}
	}

	/** Verifies that constructor and setter substring filters normalize identically.
	 * @param tree Tiny taxonomy tree
	 * @param baseline True to expect the constructor raw-case baseline failure */
	private static void checkMixedCaseContains(final TaxTree tree, final boolean baseline){
		final TaxFilter constructed=new TaxFilter(tree, -1, 0, true, true, null, null, "FIRST");
		final TaxFilter lowerConstructed=new TaxFilter(tree, -1, 0, true, true, null, null, "first");
		final TaxFilter setter=new TaxFilter(tree, true);
		setter.setContainsString("FIRST");
		final boolean constructedMatches=constructed.passesFilterByNameOnly("first leaf");
		final boolean lowerConstructedMatches=lowerConstructed.passesFilterByNameOnly("first leaf");
		final boolean setterMatches=setter.passesFilterByNameOnly("first leaf");
		check(lowerConstructedMatches, "Full constructor must retain a lower-case needle");
		check(setterMatches, "setContainsString must lowercase its needle");
		if(baseline){
			check(!constructedMatches,
					"Baseline full constructor must retain the mixed-case needle");
		}else{
			check(constructedMatches,
					"Fixed full constructor must lowercase the mixed-case needle");
		}
	}

	/** Verifies normal-return TextFile leaks from file-backed filter inputs.
	 * @param tree Tiny taxonomy tree
	 * @param baseline True to expect leaked TextFile handles
	 * @throws IOException If the test files or file-descriptor scans fail */
	private static void checkTextFileNormalReturnLeaks(final TaxTree tree,
			final boolean baseline) throws IOException{
		final Path dir=Files.createTempDirectory("taxfilter-textfile-");
		try{
			checkAddNumbersFileLeak(tree, dir.resolve("ids.txt"), baseline);
			checkReviseByBestEffortLeak(tree, dir.resolve("present.fa"), baseline);
		}finally{
			deleteTempTree(dir);
		}
	}

	/** Verifies that addNumbers closes a TextFile used for numeric ID expansion.
	 * @param tree Tiny taxonomy tree
	 * @param idFile Path for a tiny ID file
	 * @param baseline True to expect the current leak
	 * @throws IOException If the test file or file-descriptor scan fails */
	private static void checkAddNumbersFileLeak(final TaxTree tree, final Path idFile,
			final boolean baseline) throws IOException{
		Files.write(idFile, new byte[] {'4', '\n'});
		final TaxFilter filter=new TaxFilter(tree, true);
		final int before=openHandles(idFile);
		filter.addNumbers(idFile.toString(), false);
		checkTaxIDs(filter, 4);
		checkNormalReturnLeak("TaxFilter.addNumbers", idFile, before, 1, baseline);
	}

	/** Verifies that reviseByBestEffort closes its scanned FASTA TextFile.
	 * @param tree Tiny taxonomy tree
	 * @param fasta Path for a tiny FASTA file
	 * @param baseline True to expect the current leak
	 * @throws IOException If the test file or file-descriptor scan fails */
	private static void checkReviseByBestEffortLeak(final TaxTree tree,
			final Path fasta, final boolean baseline) throws IOException{
		Files.write(fasta, new byte[] {'>', 't', 'i', 'd', '|', '4', '|', '\n', 'A', '\n'});
		final TaxFilter filter=new TaxFilter(tree, TaxTree.SPECIES_E, 0, true, true,
				null, null, null);
		filter.addNumber(4, false);
		final int before=openHandles(fasta);
		filter.reviseByBestEffort(fasta.toString());
		checkTaxIDs(filter, 4);
		checkNormalReturnLeak("TaxFilter.reviseByBestEffort", fasta, before, 1,
				baseline);
	}

	/** Verifies that malformed file-expanded numeric IDs fail loudly and close.
	 * @param tree Tiny taxonomy tree
	 * @throws IOException If the test file or file-descriptor scan fails */
	private static void checkAddNumbersFileParseFailureCloses(final TaxTree tree)
			throws IOException{
		final Path dir=Files.createTempDirectory("taxfilter-malformed-");
		final Path idFile=dir.resolve("ids.txt");
		try{
			Files.write(idFile, new byte[] {'4', 'x', '\n'});
			final TaxFilter filter=new TaxFilter(tree, true);
			final int before=openHandles(idFile);
			RuntimeException failure=null;
			try{
				filter.addNumbers(idFile.toString(), false);
			}catch(RuntimeException e){
				failure=e;
			}
			check(failure instanceof NumberFormatException,
					"Malformed numeric ID file must keep NumberFormatException primary: "+
					failure);
			checkNormalReturnLeak("TaxFilter.addNumbers malformed", idFile, before,
					0, false);
		}finally{
			deleteTempTree(dir);
		}
	}

	/** Applies baseline or fixed file-handle expectations after normal return.
	 * @param label Case label
	 * @param path Scanned text file
	 * @param before Open handles before invocation
	 * @param baselineLeaks Exact number of handles expected to leak now
	 * @param baseline True to expect the current leak
	 * @throws IOException If handle counting fails */
	private static void checkNormalReturnLeak(final String label, final Path path,
			final int before, final int baselineLeaks, final boolean baseline)
			throws IOException{
		final int after=openHandles(path);
		final int leaked=after-before;
		final int expectedLeaks=(baseline ? baselineLeaks : 0);
		System.out.println(label+" textfile leaked_handles="+leaked+
				" expect="+expectedLeaks);
		check(leaked==expectedLeaks, label+" expected "+expectedLeaks+
				" leaked text-file handles but saw "+leaked+
				"; before="+before+" after="+after+" baseline="+baseline);
	}

	/** Counts file descriptors in this process that target a path.
	 * @param path Target file
	 * @return Current open file descriptor count for that exact file
	 * @throws IOException If /proc/self/fd cannot be scanned */
	private static int openHandles(final Path path) throws IOException{
		final Path target=path.toRealPath();
		int count=0;
		try(DirectoryStream<Path> stream=Files.newDirectoryStream(PROC_FD)){
			for(Path fd : stream){
				try{
					if(Files.readSymbolicLink(fd).equals(target)){count++;}
				}catch(NoSuchFileException ignored){}
			}
		}
		return count;
	}

	/** Deletes temporary test files.
	 * @param root Temporary directory root
	 * @throws IOException If cleanup fails */
	private static void deleteTempTree(final Path root) throws IOException{
		try(DirectoryStream<Path> stream=Files.newDirectoryStream(root)){
			for(Path path : stream){Files.deleteIfExists(path);}
		}
		Files.deleteIfExists(root);
	}

	/** Reports whether assertions are enabled in this JVM.
	 * @return True if assertions are enabled */
	private static boolean assertionsEnabled(){
		boolean enabled=false;
		assert(enabled=true);
		return enabled;
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Verifies an unordered set of taxonomy IDs.
	 * @param filter Filter to inspect
	 * @param expectedIDs Expected literal taxonomy IDs */
	private static void checkTaxIDs(final TaxFilter filter, final int... expectedIDs){
		final Integer[] observed=filter.taxSet();
		check(observed!=null, "Expected "+expectedIDs.length+" taxIDs but found null");
		check(observed.length==expectedIDs.length,
				"Expected "+expectedIDs.length+" taxIDs but found "+observed.length);
		for(final int expected : expectedIDs){
			boolean found=false;
			for(final Integer id : observed){
				if(id!=null && id==expected){found=true; break;}
			}
			check(found, "Missing literal taxID "+expected);
		}
	}

	/** Linux file-descriptor directory for the running fixture process. */
	private static final Path PROC_FD=new java.io.File("/proc/self/fd").toPath();
}
