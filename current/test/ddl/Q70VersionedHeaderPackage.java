package test.ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * Repackages the accepted 120-group alphabet and 2048-word reduced K6 blacklist
 * with version-one headers, preserving the mapping and loaded keys. Emits separate
 * shared/full provenance variants and a deliberately invalid version fixture.
 * Uses a fresh output directory and retains outputs for inspection. Run in a separate
 * JVM: amino-input and pigz/unpigz settings are changed and are not restored.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class Q70VersionedHeaderPackage{

	/** Repackages and checks both provenance variants, then writes summary.tsv and PASS.
	 * @param args Old alphabet, old blacklist, fresh absolute output directory,
	 * shared provenance and full provenance; the parent directory must exist */
	public static void main(String[] args) throws Exception{
		check(args.length==5, "Usage: old_alphabet old_blacklist out shared_provenance full_provenance");
		Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final Path oldAlphabet=Paths.get(args[0]), oldBlacklist=Paths.get(args[1]), out=Paths.get(args[2]);
		final String sharedProvenance=args[3], fullProvenance=args[4];
		check(out.isAbsolute() && !Files.exists(out), "Require fresh absolute output directory");
		check(sharedProvenance.length()>0 && safeHeaderValue(sharedProvenance), "Shared provenance must fit the FASTA ID header");
		check(fullProvenance.length()>0 && safeHeaderValue(fullProvenance), "Full provenance must fit the FASTA ID header");
		check(!sharedProvenance.equals(fullProvenance), "Shared and full blacklist provenances must be exact and distinct");
		Files.createDirectory(out);
		final String mapping=readMapping(oldAlphabet), aabits=mapping.substring("aabits=".length());
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6, aabits);
		final int groups=groupCount(aabits);
		check(groups==120, "Accepted R5 alphabet must retain 120 ordered-pair groups");
		final Path newAlphabet=out.resolve("alphabet.config");
		Files.write(newAlphabet, java.util.Arrays.asList("#bbgap_alphabet_config\t1", mapping), StandardCharsets.US_ASCII);
		check(readMapping(newAlphabet).equals(mapping), "Versioned alphabet must reload the exact aabits mapping");

		final ArrayList<String> lines=new ArrayList<String>(Files.readAllLines(oldBlacklist, StandardCharsets.US_ASCII));
		check(!lines.isEmpty() && lines.get(0).equals(">1"), "Accepted blacklist must start with FASTA header >1");
		final HashSet<Long> oldKeys=projectedKeys(lines, geometry);
		check(oldKeys.size()==2048, "Accepted blacklist must contain 2048 distinct reduced K6 keys");
		final Path sharedBlacklist=writeBlacklist(out.resolve("shared"), lines, oldBlacklist, oldKeys, geometry,
			sharedProvenance, aabits);
		final Path fullBlacklist=writeBlacklist(out.resolve("full"), lines, oldBlacklist, oldKeys, geometry,
			fullProvenance, aabits);
		Files.write(out.resolve("summary.tsv"), java.util.Arrays.asList("file\tvalue",
			"alphabet\t"+newAlphabet,
			"shared_blacklist\t"+sharedBlacklist,
			"full_blacklist\t"+fullBlacklist,
			"words\t2048",
			"alphabet_groups\t"+groups,
			"shared_provenance\t"+sharedProvenance,
			"full_provenance\t"+fullProvenance), StandardCharsets.US_ASCII);
		Files.write(out.resolve("PASS"), java.util.Arrays.asList(
			"Q70_VERSIONED_HEADER_PACKAGE_PASS words=2048 alphabet_groups=120"), StandardCharsets.US_ASCII);
		System.out.println("Q70_VERSIONED_HEADER_PACKAGE_PASS words=2048 alphabet_groups=120");
	}

	/*--------------------------------------------------------------*/
	/*----------------     Resource conversion      ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies the original lines, replaces only their first header and verifies exact loaded-key equality.
	 * Retains an unknown-version rejection fixture beside the valid output; returns the valid path.
	 * The supplied list and key set are borrowed and unchanged. */
	private static Path writeBlacklist(Path out, ArrayList<String> original, Path oldBlacklist, HashSet<Long> oldKeys,
			ProteinKmerGeometry geometry, String provenance, String aabits) throws Exception{
		Files.createDirectory(out);
		final ArrayList<String> lines=new ArrayList<String>(original);
		final String header="bbgap_blacklist_v1;k=6;words=2048;mode=reduced_literal;provenance="+headerValue(provenance)+";aabits="+aabits;
		lines.set(0, ">"+header);
		final Path newBlacklist=out.resolve("reduced_shadow_k6_2048.faa");
		Files.write(newBlacklist, lines, StandardCharsets.US_ASCII);
		check(projectedKeys(Files.readAllLines(newBlacklist, StandardCharsets.US_ASCII), geometry).equals(oldKeys),
			"Versioned blacklist FASTA must preserve the exact reduced K6 key set");
		final DDLPentamerBlacklist oldFilter=DDLPentamerBlacklist.loadReducedStructural(oldBlacklist, "unmeasured",
			2048, geometry, provenance);
		final DDLPentamerBlacklist newFilter=DDLPentamerBlacklist.loadReducedStructural(newBlacklist, "unmeasured",
			2048, geometry, provenance);
		check(oldFilter.size()==2048 && newFilter.size()==2048, "Both blacklist readers must load exactly 2048 words");
		for(long key:oldKeys){
			check(oldFilter.contains(key) && newFilter.contains(key),
				"Versioned blacklist load must preserve every accepted reduced K6 key");
		}

		final ArrayList<String> bad=new ArrayList<String>(lines);
		bad.set(0, ">"+header.replace("bbgap_blacklist_v1", "bbgap_blacklist_v2"));
		final Path unknown=out.resolve("reduced_shadow_k6_2048.unknown_version.faa");
		Files.write(unknown, bad, StandardCharsets.US_ASCII);
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(unknown, "unmeasured", 2048, geometry, provenance));
		return newBlacklist;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Literal input helpers     ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads one aabits line, allowing empty lines and an optional unique first-line v1 header.
	 * Other comment headers, multiple mappings or a missing mapping fail. */
	private static String readMapping(Path path) throws Exception{
		String mapping=null;
		int headers=0;
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		for(int i=0; i<lines.size(); i++){
			final String line=lines.get(i);
			if(line.equals("#bbgap_alphabet_config\t1")){
				check(i==0 && mapping==null && ++headers==1, "Alphabet version header must be first and unique");
			}else if(line.startsWith("#")){
				throw new IllegalArgumentException("Unknown alphabet header: "+line);
			}else if(line.length()>0){
				check(mapping==null && line.startsWith("aabits="), "Alphabet config must contain one aabits line");
				mapping=line;
			}
		}
		check(mapping!=null, "Alphabet config is missing aabits");
		return mapping;
	}

	/** Requires the dimer:1,3,5: prefix and counts slash-separated groups, including empty fields. */
	private static int groupCount(String aabits){
		final String prefix="dimer:1,3,5:";
		check(aabits.startsWith(prefix), "Q70 blacklist uses the explicit dimer K6 alphabet");
		return aabits.substring(prefix.length()).split("/", -1).length;
	}

	/** Rejects Character.isWhitespace characters; callers separately require a nonempty value. */
	private static boolean safeHeaderValue(String value){
		for(int i=0; i<value.length(); i++){
			if(Character.isWhitespace(value.charAt(i))){return false;}
		}
		return true;
	}

	/** Percent-encodes percent and semicolon delimiters, retaining other characters unchanged. */
	private static String headerValue(String value){
		final StringBuilder encoded=new StringBuilder(value.length());
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c=='%' || c==';'){
				encoded.append('%').append(Character.forDigit((c>>>4)&15, 16)).append(Character.forDigit(c&15, 16));
			}else{encoded.append(c);}
		}
		return encoded.toString();
	}

	/** Joins nonempty sequence lines between FASTA headers and returns distinct projected K6 keys.
	 * The caller supplies the native fixture lines; duplicate projected keys are rejected. */
	private static HashSet<Long> projectedKeys(List<String> lines, ProteinKmerGeometry geometry){
		final HashSet<Long> keys=new HashSet<Long>();
		final StringBuilder bases=new StringBuilder();
		boolean seen=false;
		for(String line:lines){
			if(line.startsWith(">")){
				if(seen){add(keys, bases, geometry);}
				seen=true;
				bases.setLength(0);
			}else if(line.length()>0){bases.append(line.trim());}
		}
		if(seen){add(keys, bases, geometry);}
		return keys;
	}

	/** Requires six canonical uppercase residues, packs their five-bit codes and inserts one reduced key. */
	private static void add(HashSet<Long> keys, StringBuilder bases, ProteinKmerGeometry geometry){
		check(bases.length()==6, "Blacklist records must be K6");
		long raw=0;
		final String aa="ARNDCEQGHILKMFPSTWYV";
		for(int i=0; i<bases.length(); i++){
			final int value=aa.indexOf(bases.charAt(i));
			check(value>=0, "Blacklist records must use canonical uppercase amino acids");
			raw=(raw<<5)|value;
		}
		check(keys.add(geometry.word(raw)), "Versioned packaging must preserve a distinct reduced key set");
	}

	/** Requires IllegalArgumentException; checked exceptions do not count as successful rejection. */
	private static void expectFailure(Action action) throws Exception{
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Expected malformed versioned blacklist to fail");
	}

	/** Throws IllegalArgumentException on a failed fixture invariant, independent of -ea. */
	private static void check(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** One checked fixture operation, permitting I/O exceptions to propagate. */
	private interface Action{
		/** Executes the malformed-input operation. */
		void run() throws Exception;
	}
}
