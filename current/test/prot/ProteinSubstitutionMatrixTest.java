package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;

import dna.AminoAcid;
import prot.Blosum62;
import prot.FamilyAlignmentMetric;
import prot.FamilyAlignmentMetric.Mode;
import prot.FamilyAlignmentMetric.Observation;
import prot.HbmPositionModel;
import prot.ProteinSequence;
import prot.ProteinSubstitutionMatrix;
import prot.ProteinSubstitutionMatrixIO;
import prot.cluster.TailAlignmentGate;

/** Literal score/trace witnesses and independent default-versus-scaled controls.
 * Retains generated evidence in the supplied existing directory and rejects overwrites.
 * Uses synthetic inputs plus the installed accepted default matrix, not a biological panel.
 * Run in an isolated JVM: evidence bookkeeping is static. The private resource override
 * is restored after its rejection checks; worker-isolation checks own and join their threads.
 * @author Brian, Keqing, Fischl
 */
public final class ProteinSubstitutionMatrixTest{

	/** Runs native/scaled/runtime/interchange checks, then writes PASS and an explicit file roster.
	 * @param args One existing fixture directory; all output names used here must be absent */
	public static void main(String[] args)throws Exception{
		check(args.length==1 && Files.isDirectory(Paths.get(args[0])), "Supply an existing native fixture directory");
		final Path root=Paths.get(args[0]).toAbsolutePath().normalize();
		evidenceRoot=root;
		evidence.clear();
		final ProteinSubstitutionMatrix baseline=ProteinSubstitutionMatrix.BLOSUM62;
		final int[][] table=scaledScores(4);
		final ProteinSubstitutionMatrix scaled=ProteinSubstitutionMatrix.fromScores("custom40", 4, 16, table);
		table[0][0]=99;
		check(scaled.score((byte)0, (byte)0)==4*Blosum62.score((byte)0, (byte)0), "Matrix must detach caller-owned rows");
		check(scaled.score(Blosum62.X_CODE, Blosum62.X_CODE)==-4, "X retains scaled legacy policy");
		check(scaled.same(ProteinSubstitutionMatrix.parseBinding(scaled.binding())),
			"Portable binding must preserve every score and unit");
		check(!baseline.same(scaled) && !baseline.sha80.equals(scaled.sha80),
			"Calibration identity must distinguish scaled matrices");
		final int[][] altered=scaledScores(4);
		altered[0][1]++;
		altered[1][0]++;
		final ProteinSubstitutionMatrix sameName=ProteinSubstitutionMatrix.fromScores("custom40", 4, 16, altered);
		check(!scaled.same(sameName) && !scaled.sha80.equals(sameName.sha80) &&
			!scaled.same(ProteinSubstitutionMatrix.parseBinding(sameName.binding())) &&
			sameName.same(ProteinSubstitutionMatrix.parseBinding(sameName.binding())),
			"Same name and units must not hide different score bytes, including after portable reload");
		check(ProteinSubstitutionMatrix.load("blosum62", null)==baseline,
			"Explicit legacy selector must retain the built-in object");
		final ProteinSubstitutionMatrix adopted=ProteinSubstitutionMatrix.defaultMatrix();
		check(adopted.same(ProteinSubstitutionMatrix.load("POSSUM55", null)) &&
			adopted.same(ProteinSubstitutionMatrix.read(ProteinSubstitutionMatrix.defaultResourcePath())) &&
			adopted.name.equals("POSSUM55") && adopted.sha80.equals("21a59255f7b0dfa77945") &&
			adopted.scale==4 && adopted.linearGap==16,
			"New default must load the accepted portable POSSUM55 scores and units");
		final String defaultPath=ProteinSubstitutionMatrix.defaultResourcePath();
		check(Paths.get(defaultPath).getFileName().toString().equals("POSSUM55.tsv"),
			"Installed release lookup must use the uppercase filename on case-sensitive filesystems");
		check(prot.DigestSuffix.file(defaultPath).equals("474b6013e399137f0fa1"),
			"Installed default bytes must match the accepted resource");
		check(adopted.score((byte)aa('A'), (byte)aa('A'))==17 &&
			adopted.score((byte)aa('A'), (byte)aa('C'))==-1,
			"Accepted POSSUM55 witness scores distinguish it from the legacy custom40 matrix");
		final String originalOverride=System.getProperty("bbgap.custom40.resource");
		try{
			for(String invalid:new String[]{"", "relative.tsv", root.toAbsolutePath().toString(),
				root.toAbsolutePath().resolve("..").resolve("unresolved.tsv").toString()}){
				System.setProperty("bbgap.custom40.resource", invalid);
				check(defaultPath.equals(ProteinSubstitutionMatrix.defaultResourcePath()),
					"A legacy custom40 override must not redirect the POSSUM55 default");
				expectFailure(()->ProteinSubstitutionMatrix.load("custom40", null),
					"Invalid explicit legacy resource override must fail without poisoning later loads");
			}
		}finally{
			if(originalOverride==null){
				System.clearProperty("bbgap.custom40.resource");
			}else{
				System.setProperty("bbgap.custom40.resource", originalOverride);
			}
		}
		final ProteinSubstitutionMatrix legacyDefault=ProteinSubstitutionMatrix.load("custom40", null);
		check(legacyDefault.sha80.equals("f3ddb77a31c0b04e56dc") && !legacyDefault.same(adopted) &&
			legacyDefault.score((byte)aa('A'), (byte)aa('A'))==15 &&
			legacyDefault.score((byte)aa('A'), (byte)aa('C'))==2,
			"Explicit custom40 selection retains its distinct historical scores after the default migration");
		final java.util.HashMap<String, String> historical=new java.util.HashMap<String, String>();
		check(prot.FamilyScoringPackage.resolveMatrix(historical)==baseline,
			"Historical absent-selector package must remain BLOSUM after default adoption");
		historical.put("matrix", "custom40");
		historical.put("matrixbinding", scaled.binding());
		check(prot.FamilyScoringPackage.resolveMatrix(historical).same(scaled) && !scaled.same(adopted),
			"An explicit historical binding retains its actual scores even when the name matches today's default");
		expectFailure(()->ProteinSubstitutionMatrix.load("unknown", null), "Unknown selector must fail");
		expectFailure(()->ProteinSubstitutionMatrix.load("custom60", null), "Nondefault custom resource remains required");
		expectFailure(()->ProteinSubstitutionMatrix.parseBinding("custom40:4:16:1"), "Truncated binding must fail");
		expectFailure(()->ProteinSubstitutionMatrix.fromScores("blosum62", 4, 16, scaledScores(4)),
			"Reserved name must not shadow BLOSUM62");
		final int[][] asymmetric=scaledScores(4);
		asymmetric[0][1]++;
		expectFailure(()->ProteinSubstitutionMatrix.fromScores("custom40", 4, 16, asymmetric), "Asymmetric scores must fail");
		final int[][] diagonal=scaledScores(4);
		diagonal[0][0]=0;
		expectFailure(()->ProteinSubstitutionMatrix.fromScores("custom40", 4, 16, diagonal), "Zero self diagonal must fail");
		final int[][] oversized=scaledScores(4);
		oversized[0][0]=101;
		expectFailure(()->ProteinSubstitutionMatrix.fromScores("custom40", 4, 16, oversized),
			"Scores beyond the fixed DP bound must fail");

		final String resource=resource(scaled);
		final Path valid=root.resolve("custom40.tsv");
		write(valid, resource);
		check(scaled.same(ProteinSubstitutionMatrix.load("custom40", valid.toString())),
			"TSV letter mapping must round-trip all400 cells");
		expectFailure(()->ProteinSubstitutionMatrix.load("custom60", valid.toString()),
			"Selector must agree with resource identity");
		final Path duplicate=root.resolve("duplicate.tsv");
		write(duplicate, resource+"score\tA\tA\t16\n");
		expectFailure(()->ProteinSubstitutionMatrix.read(duplicate.toString()), "Repeated score row must fail");
		final Path missing=root.resolve("missing.tsv");
		write(missing, resource.substring(0, resource.lastIndexOf("score\t")));
		expectFailure(()->ProteinSubstitutionMatrix.read(missing.toString()), "Missing score row must fail");
		final Path malformed=root.resolve("malformed.tsv");
		write(malformed, resource.replace("score\tA\tA\t16", "score\tA\tA\tNaN"));
		expectFailure(()->ProteinSubstitutionMatrix.read(malformed.toString()), "Noninteger score must fail");
		releaseName(root, scaled);

		for(Mode mode:new Mode[]{Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, Mode.PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY}){
			final FamilyAlignmentMetric ordinary=new FamilyAlignmentMetric(mode, 3);
			final FamilyAlignmentMetric explicit=new FamilyAlignmentMetric(mode, 3, baseline);
			final FamilyAlignmentMetric precise=new FamilyAlignmentMetric(mode, 3, scaled);
			reuse(precise);
			check(ordinary.same(explicit) && !ordinary.same(precise), "Metric equality includes immutable matrix identity");
			for(String query:new String[]{"AAA", "RAAA", "AAAR", "AARAA", "ACDEFGHIK", "RRR"}){
				for(String reference:new String[]{"AAA", "AAAA", "ACDEFGHIK", "RRR"}){
					final ProteinSequence a=sequence(query), b=sequence(reference);
					final int oldSelf=baseline.selfScore(a.enc);
					final Observation x=new Observation(), y=new Observation(), z=new Observation();
					ordinary.newWorker().measure(a, b, oldSelf, x);
					explicit.newWorker().measure(a, b, oldSelf, y);
					precise.newWorker().measure(a, b, oldSelf, z);
					compare(x, y, 1);
					compare(x, z, 4);
				}
			}
		}
		final byte[] gapQuery=Blosum62.encode("AARAA".getBytes(StandardCharsets.US_ASCII), "fixture");
		final byte[] gapReference=Blosum62.encode("AAAA".getBytes(StandardCharsets.US_ASCII), "fixture");
		final HbmPositionModel ordinary=HbmPositionModel.fromConsensus(gapReference);
		final HbmPositionModel precise=HbmPositionModel.fromConsensus(gapReference, scaled);
		final HbmPositionModel.Result ox=ordinary.align(gapQuery, true), sx=precise.align(gapQuery, true);
		check(ox.score==12*64 && sx.score==48*64 && Arrays.equals(ox.path, sx.path),
			"Four A pairs minus one linear gap must give12/48 and the same glocal trace");
		check(precise.scorePath(gapQuery, sx.start, sx.path)==sx.score, "Custom trace rescoring must use the same gap units");
		check(precise.scoreTail(gapQuery, 5)==4*ordinary.scoreTail(gapQuery, 5),
			"Rolling tail scorer must share the selected gap");

		final int[][] changedScores=scaledScores(4);
		final int a=AminoAcid.acidToNumberExtended['A'], r=AminoAcid.acidToNumberExtended['R'];
		changedScores[a][r]=changedScores[r][a]=16;
		final ProteinSubstitutionMatrix changed=ProteinSubstitutionMatrix.fromScores("custom60", 4, 16, changedScores);
		final FamilyAlignmentMetric selected=new FamilyAlignmentMetric(Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, 3, changed);
		final ProteinSequence query=sequence("RAAA"), reference=sequence("AAAA");
		final Observation actual=new Observation();
		selected.newWorker().measure(query, reference, 17, actual);
		check(actual.measured() && actual.rawScore==64 && actual.denominatorSelf==68 && actual.fullQuerySelf==68 &&
			actual.queryStart==0 && actual.referenceStart==0 && actual.columns==4 && actual.identities==3 &&
			actual.identityPct()==75 && actual.blosum()==64.0/68,
			"Positive custom R/A score selects all four pairs and the custom diagonal denominator");
		final TailAlignmentGate gate=new TailAlignmentGate(TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT,
			TailAlignmentGate.Denominator.ALIGNED_QUERY, 4, 3, 0, 0, null, changed);
		final HbmPositionModel.ConsensusScorer wrong=new HbmPositionModel.ConsensusScorer().resetConsensus(reference.enc);
		expectFailure(()->gate.decide(query.enc, reference.enc, wrong), "Trace scorer and gate matrix mismatch must fail");
		expectFailure(()->new FamilyAlignmentMetric(Mode.GLOCAL_FULL_QUERY, 0, changed),
			"Unsupported legacy custom package must fail explicitly");
		parallel(query, reference, selected, new FamilyAlignmentMetric(Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, 3));
		interchange(root);
		finish("PROTEIN_SUBSTITUTION_MATRIX_TEST_PASS loader=true malformed=true default_parity=true scale_parity=true custom_trace=true worker_isolation=true");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Text Interchange       ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses a literal nonuniform distribution and hand-calculated score oracles.
	 * Native/runtime identity is checked separately from exported ambiguity values.
	 */
	private static void interchange(Path root)throws Exception{
		final int[][] scores=new int[20][20];
		for(int i=0; i<20; i++){
			Arrays.fill(scores[i], -4);
			scores[i][i]=8;
		}
		pair(scores, 'A', 'N', 2);
		pair(scores, 'A', 'E', 2);
		pair(scores, 'A', 'Q', 4);
		final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.fromScores("POSSUM55", 4, 16, scores);
		final double[] background=new double[20];
		background[aa('A')]=.25;
		background[aa('D')]=.375;
		background[aa('N')]=.125;
		background[aa('E')]=.125;
		background[aa('Q')]=.125;
		final double[] original=background.clone();
		final String text=ProteinSubstitutionMatrixIO.ncbiText(matrix, background, "literal fixture").toString();
		check(Arrays.equals(background, original), "Export must preserve borrowed background probabilities");
		check(ncbiCell(text, 'B', 'A')==-3, "D/N=.75/.25 gives -2.5, rounded away from zero to -3");
		check(ncbiCell(text, 'Z', 'A')==3, "Equal E/Q weights average A scores2 and4 to3");
		check(ncbiCell(text, 'X', 'A')==2, "Background expectation against A is1.5, rounded to2");
		check(ncbiCell(text, 'B', 'B')==4, "Joint D/N expectation is3.5, rounded once to4");
		check(ncbiCell(text, 'Z', 'Z')==2 && ncbiCell(text, 'B', 'Z')==-4,
			"Joint ambiguity expectations preserve literal Z/Z=2 and B/Z=-4");
		check(ncbiCell(text, 'B', 'X')==0 && ncbiCell(text, 'Z', 'X')==-1 && ncbiCell(text, 'X', 'X')==0,
			"Direct unrounded joint expectations .125, -.75 and .25 round to0,-1,0");
		for(char a:"ARNDCQEGHILKMFPSTWYVBZX*".toCharArray()){
			for(char b:"ARNDCQEGHILKMFPSTWYVBZX*".toCharArray()){
				check(ncbiCell(text, a, b)==ncbiCell(text, b, a), "Export must be exactly symmetric");
			}
			check(ncbiCell(text, '*', a)==(a=='*' ? 1 : -4), "Star uses the literal matrix minimum or self+1");
		}
		check(text.contains("16 score units") && text.contains("Statistical E-value parameters are not provided"),
			"Generic export comments state actual gap units and statistical limits");
		check(!text.contains("Identity band:") && !text.contains("8*log2"),
			"A POSSUM55 name alone must not turn a synthetic matrix into a claimed fitted resource");
		final Path exported=root.resolve("ncbi-disguised.tsv");
		write(exported, text);
		final ProteinSubstitutionMatrix restored=ProteinSubstitutionMatrix.read(exported.toString());
		check(matrix.same(restored), "NCBI round-trip preserves name, scale, gap and every canonical score");
		check(restored.score(Blosum62.X_CODE, Blosum62.X_CODE)==-4,
			"Exported X/X=0 must not replace the existing runtime X/X=-4");
		for(byte a=0; a<22; a++){
			if(a==20){continue;}
			for(byte b=0; b<22; b++){
				if(b!=20){check(restored.score(a, b)==matrix.score(a, b), "Every runtime score remains unchanged");}
			}
		}
		final byte[] encoded=new ProteinSequence("ambiguous", "ABZJX*").enc;
		check(encoded.length==5 && encoded[1]==Blosum62.X_CODE && encoded[2]==Blosum62.X_CODE
			&& encoded[3]==Blosum62.X_CODE && encoded[4]==Blosum62.X_CODE,
			"Interchange preserves ambiguity folding and one-terminal-stop handling");
		check(matrix.selfScore(encoded)==restored.selfScore(encoded), "Round-trip self scores remain identical");
		final byte[] reference=sequence("AAAADN").enc, query=sequence("AADNX").enc;
		final HbmPositionModel.Result before=HbmPositionModel.fromConsensus(reference, matrix).align(query, true);
		final HbmPositionModel.Result after=HbmPositionModel.fromConsensus(reference, restored).align(query, true);
		check(before.score==after.score && before.start==after.start && Arrays.equals(before.path, after.path),
			"Round-trip preserves consensus score, start and complete traceback");

		final String nativeText=ProteinSubstitutionMatrixIO.nativeText(matrix).toString();
		check(nativeText.equals(resource(matrix)), "Native writer exactly matches the pre-existing literal serializer");
		final Path nativeFile=root.resolve("native-disguised.mat");
		write(nativeFile, nativeText);
		check(matrix.same(ProteinSubstitutionMatrix.read(nativeFile.toString())), "Native detection ignores filename suffix");
		final Path plain=root.resolve("plain-reordered.ncbi");
		final String plainText=ncbiFixture(matrix, false);
		write(plain, plainText);
		check(matrix.same(ProteinSubstitutionMatrixIO.read(plain.toString(), "POSSUM55", 4, 16)),
			"Reordered rows/columns, whitespace and ignored integer J/ambiguity cells preserve canonical scores");
		expectFailure(()->ProteinSubstitutionMatrix.read(plain.toString()), "explicit name, scale and gap",
			"Plain input cannot infer units or name from its file");
		expectFailure(()->ProteinSubstitutionMatrixIO.read(exported.toString(), "custom40", 4, 16),
			"Conflicting matrix name", "Explicit name cannot override embedded metadata");
		expectFailure(()->ProteinSubstitutionMatrixIO.read(exported.toString(), "POSSUM55", 8, 16),
			"Conflicting matrix scale", "Eighth-bit units do not mean native scale8");
		expectFailure(()->ProteinSubstitutionMatrixIO.read(nativeFile.toString(), "POSSUM55", 4, 4),
			"Conflicting matrix gap", "Explicit metadata must agree with native metadata too");
		malformedNcbi(root, "bad-extra", plainText.replaceFirst("999", "nope"), "Noninteger");
		malformedNcbi(root, "overflow-extra", plainText.replaceFirst("999", "2147483648"), "overflow");
		malformedNcbi(root, "unknown-symbol", plainText.replace('J', 'U'), "Unknown NCBI residue");
		malformedNcbi(root, "duplicate-column", plainText.replaceFirst("J ", "B "), "Duplicate NCBI column");
		malformedNcbi(root, "asymmetric", ncbiFixture(matrix, true), "symmetric");
		final int last=plainText.lastIndexOf('\n', plainText.length()-2);
		malformedNcbi(root, "missing-row", plainText.substring(0, last+1), "Incomplete NCBI");
		malformedNcbi(root, "duplicate-row", plainText+plainText.substring(last+1), "duplicate NCBI row");
		final Path badNative=root.resolve("bad-native-version.mat");
		write(badNative, nativeText.replace("#bbgap_substitution_matrix\t1", "#bbgap_substitution_matrix\t2"));
		expectFailure(()->ProteinSubstitutionMatrix.read(badNative.toString()), "Unexpected substitution matrix schema",
			"Unknown native version must not fall through to NCBI parsing");
		final Path displaced=root.resolve("displaced-native.mat");
		write(displaced, "\n"+nativeText);
		expectFailure(()->ProteinSubstitutionMatrix.read(displaced.toString()), "exact first line",
			"Native header remains strictly first, even though NCBI supports blank lines");
		final double[] nonfinite=background.clone();
		nonfinite[0]=Double.NaN;
		expectFailure(()->ProteinSubstitutionMatrixIO.ncbiText(matrix, nonfinite, "bad"), "finite",
			"Nonfinite background must fail before output");
		final double[] noDn=new double[20];
		noDn[aa('E')]=.5;
		noDn[aa('Q')]=.5;
		expectFailure(()->ProteinSubstitutionMatrixIO.ncbiText(matrix, noDn, "bad"), "Undefined",
			"Undefined D/N weights cannot be replaced by uniform guesses");
		final double[] unnormalized=background.clone();
		unnormalized[aa('A')]=.5;
		expectFailure(()->ProteinSubstitutionMatrixIO.ncbiText(matrix, unnormalized, "bad"), "sum to one",
			"Background mass must be explicit and normalized");
		derivedConversion(root, nativeFile);
		write(root.resolve("NCBI_PASS"), "NCBI_MATRIX_INTERCHANGE_PASS canonical=400 literal_extras=true native_parity=true runtime_parity=true malformed=true\n");
	}

	/** A synthetic directed table has literal diagonal27/off-diagonal-8 scores.
	 * Diagonal probability21/800 and off-diagonal1/800 give uniform endpoints1/20;
	 * this fixture verifies pairing and metadata, not biological provenance.
	 */
	private static void derivedConversion(Path root, Path mismatching)throws Exception{
		final int[][] scores=new int[20][20];
		for(int i=0; i<20; i++){
			Arrays.fill(scores[i], -8);
			scores[i][i]=27;
		}
		final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.fromScores("POSSUM55", 4, 16, scores);
		final Path source=root.resolve("derived-literal.tsv"), background=root.resolve("literal-probabilities.tsv");
		write(source, resource(matrix));
		final String labels="ARNDCEQGHILKMFPSTWYV";
		final StringBuilder state=new StringBuilder("#bbgap_substitution_probabilities\t1\nkind\tmonomer\n");
		state.append("band\tquery_symbol\tconsensus_symbol\tprobability\tquery_marginal\tconsensus_marginal\n");
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){
				state.append("25to55\t").append(labels.charAt(a)).append('\t').append(labels.charAt(b))
					.append('\t').append(a==b ? "0.02625" : "0.00125").append("\t0.05\t0.05\n");
			}
		}
		write(background, state.toString());
		final Path exported=root.resolve("derived-literal.mat"), roundtrip=root.resolve("derived-roundtrip.tsv");
		ProteinSubstitutionMatrixIO.main(new String[]{"in="+source, "out="+exported, "outformat=ncbi",
			"background="+background, "band=25to55"});
		evidence.add(root.relativize(exported).toString());
		final String description=new String(Files.readAllBytes(exported), StandardCharsets.US_ASCII);
		check(description.contains("Prokaryotic Ortholog Sequence SUbstitution Matrix")
			&& description.contains("25-55%") && description.contains("8*log2") && description.contains("eighth-bit")
			&& description.contains("Supplied25to55 probability state reproduces every canonical score"),
			"Verified matrix/state pairing permits the explicit POSSUM55 derivation description");
		ProteinSubstitutionMatrixIO.main(new String[]{"in="+exported, "out="+roundtrip, "outformat=native"});
		evidence.add(root.relativize(roundtrip).toString());
		check(Arrays.equals(Files.readAllBytes(source), Files.readAllBytes(roundtrip)),
			"Actual conversion CLI preserves every native byte through its own NCBI output");
		final Path rejected=root.resolve("mismatched-background-output.mat");
		expectFailure(()->ProteinSubstitutionMatrixIO.main(new String[]{"in="+mismatching, "out="+rejected,
			"background="+background, "band=25to55"}), "matrix/background mismatch",
			"An arbitrary POSSUM55 name must not claim derivation from an unrelated probability table");
		check(!Files.exists(rejected), "Mismatched matrix/background must fail before creating output");
		expectFailure(()->ProteinSubstitutionMatrixIO.main(new String[]{"in="+source, "out="+source,
			"outformat=native", "ow=t"}), "differ", "Overwrite permission does not permit destroying input");
		check(Arrays.equals(Files.readAllBytes(source), Files.readAllBytes(roundtrip)),
			"Rejected alias conversion preserves the source file");
	}

	/** Assigns a literal symmetric fixture pair in encoded order. */
	private static void pair(int[][] scores, char a, char b, int value){
		scores[aa(a)][aa(b)]=scores[aa(b)][aa(a)]=value;
	}

	/** Maps only canonical fixture letters; the production parser is not used as an oracle. */
	private static int aa(char c){return AminoAcid.acidToNumberExtended[c];}

	/** Independently reads one exported cell through its labels for literal-value checks. */
	private static int ncbiCell(String text, char row, char column){
		String[] labels=null;
		for(String line:text.split("\n")){
			line=line.trim();
			if(line.isEmpty() || line.startsWith("#")){continue;}
			final String[] fields=line.split("\\s+");
			if(labels==null){labels=fields;}else if(fields[0].equals(Character.toString(row))){
				for(int i=0; i<labels.length; i++){
					if(labels[i].equals(Character.toString(column))){return Integer.parseInt(fields[i+1]);}
				}
			}
		}
		throw new AssertionError("Missing exported literal cell "+row+"/"+column);
	}

	/** Emits unrelated input-axis orders and conspicuous ignored extra scores.
	 * The asymmetry switch creates one invalid canonical cell, independently of the reader.
	 */
	private static String ncbiFixture(ProteinSubstitutionMatrix matrix, boolean asymmetric){
		final String labels="J*XZBVYWTSPFMKLIHGEQCDNRA";
		final StringBuilder out=new StringBuilder("# Plain external fixture\n \t\n  ");
		for(int i=0; i<labels.length(); i++){out.append(labels.charAt(i)).append(" \t");}
		out.append('\n');
		for(int a=labels.length()-1; a>=0; a--){
			final char row=labels.charAt(a);
			out.append(row);
			for(int b=0; b<labels.length(); b++){
				final char column=labels.charAt(b);
				final boolean canonical="BZX*J".indexOf(row)<0 && "BZX*J".indexOf(column)<0;
				int value=canonical ? matrix.score((byte)aa(row), (byte)aa(column)) : 999;
				if(asymmetric && row=='A' && column=='R'){value++;}
				out.append(' ').append(value);
			}
			out.append('\n');
		}
		return out.toString();
	}

	/** Retains an invalid input and requires rejection at the intended diagnostic. */
	private static void malformedNcbi(Path root, String label, String text, String diagnostic)throws Exception{
		final Path file=root.resolve("ncbi-"+label+".mat");
		write(file, text);
		expectFailure(()->ProteinSubstitutionMatrixIO.read(file.toString(), "POSSUM55", 4, 16),
			diagnostic, "Malformed NCBI fixture "+label+" must fail");
	}

	/*--------------------------------------------------------------*/
	/*----------------   Runtime identity checks    ----------------*/
	/*--------------------------------------------------------------*/

	/** Naming changes identity, not supplied scores; explicit metadata remains strict. */
	private static void releaseName(Path root, ProteinSubstitutionMatrix legacy)throws Exception{
		final ProteinSubstitutionMatrix named=ProteinSubstitutionMatrix.fromScores("POSSUM55", 4, 16, scaledScores(4));
		check(ProteinSubstitutionMatrix.supportedSelector("POSSUM55"), "The release selector must be admitted explicitly");
		check(!ProteinSubstitutionMatrix.supportedSelector("possum55"),
			"The release has one case-sensitive selector, without a lowercase alias");
		expectFailure(()->ProteinSubstitutionMatrix.fromScores("possum55", 4, 16, scaledScores(4)),
			"exactly POSSUM55", "A lowercase resource header must not create a second release name");
		expectFailure(()->ProteinSubstitutionMatrix.parseBinding("possum55"+named.binding().substring(8)),
			"exactly POSSUM55", "A lowercase embedded name must not bypass the selector rule");
		check(!named.same(legacy) && !named.sha80.equals(legacy.sha80),
			"Renaming must not erase name-sensitive calibration identity");
		check(named.scale==4 && named.linearGap==16, "Release naming must preserve score and gap units");
		for(int a=0; a<22; a++){
			if(a==20){continue;}
			for(int b=0; b<22; b++){
				if(b==20){continue;}
				check(named.score((byte)a, (byte)b)==4*Blosum62.score((byte)a, (byte)b),
					"Release naming must preserve every canonical and ambiguity score: "+a+","+b);
			}
		}
		check(ProteinSubstitutionMatrix.load("POSSUM55", null).same(ProteinSubstitutionMatrix.defaultMatrix()) &&
			!named.same(ProteinSubstitutionMatrix.defaultMatrix()),
			"The implicit release selector loads the accepted resource, not a same-named synthetic matrix");
		final String plain=resource(named);
		final String description="description\tProkaryotic Ortholog Sequence SUbstitution Matrix\n";
		final String described=plain.replace("name\tPOSSUM55\n", "name\tPOSSUM55\n"+description);
		final Path bare=root.resolve("POSSUM55-bare.tsv"), file=root.resolve("POSSUM55.tsv");
		write(bare, plain);
		write(file, described);
		check(named.same(ProteinSubstitutionMatrix.load("POSSUM55", bare.toString())),
			"Description must remain optional for existing native files");
		check(named.same(ProteinSubstitutionMatrix.load("POSSUM55", file.toString())),
			"Description must not alter the complete score binding");
		check(named.same(ProteinSubstitutionMatrix.parseBinding(named.binding())),
			"Embedded release binding must retain its name and all scores");
		final java.util.HashMap<String, String> options=new java.util.HashMap<String, String>();
		options.put("matrix", "POSSUM55");
		options.put("matrixfile", file.toString());
		check(named.same(prot.FamilyScoringPackage.resolveMatrix(options)),
			"Package file resolution must admit the release selector");
		options.remove("matrixfile");
		options.put("matrixbinding", named.binding());
		check(named.same(prot.FamilyScoringPackage.resolveMatrix(options)),
			"Package embedded resolution must admit the same release binding");
		expectFailure(()->ProteinSubstitutionMatrix.load("custom40", file.toString()), "selector/file name mismatch",
			"Renamed resource must not masquerade as custom40");
		final Path duplicate=root.resolve("POSSUM55-description-duplicate.tsv");
		write(duplicate, described+description);
		expectFailure(()->ProteinSubstitutionMatrix.read(duplicate.toString()), "Repeated matrix description",
			"Duplicate descriptions must fail at their own guard");
		final Path wrong=root.resolve("POSSUM55-description-wrong.tsv");
		write(wrong, described.replace(description, "description\tDifferent matrix\n"));
		expectFailure(()->ProteinSubstitutionMatrix.read(wrong.toString()), "Unsupported matrix description",
			"Unreviewed description must not be silently ignored");
		final Path misplaced=root.resolve("custom40-description.tsv");
		write(misplaced, described.replace("name\tPOSSUM55\n", "name\tcustom40\n"));
		expectFailure(()->ProteinSubstitutionMatrix.read(misplaced.toString()), "description is reserved for POSSUM55",
			"Release description must not relabel a historical matrix");
		final Path unknown=root.resolve("POSSUM55-unknown-metadata.tsv");
		write(unknown, plain+"unrecognized\tvalue\n");
		expectFailure(()->ProteinSubstitutionMatrix.read(unknown.toString()), "Unknown matrix metadata",
			"New description support must not admit arbitrary metadata");
	}

	/** Reuses custom capacity across growth, shrinkage and a span rejection, comparing fresh workers. */
	private static void reuse(FamilyAlignmentMetric metric){
		final FamilyAlignmentMetric.Worker reused=metric.newWorker();
		final String[] queries={"AAA", "ACDEFGHIKLMNPQRSTVWY", "RAAA", "AA", "AAAR"};
		final String[] references={"AAA", "AACDEFGHIKLMNPQRSTVWY", "AAAA", "AAA", "AAA"};
		for(int i=0; i<queries.length; i++){
			final ProteinSequence q=sequence(queries[i]), r=sequence(references[i]);
			final int self=ProteinSubstitutionMatrix.BLOSUM62.selfScore(q.enc);
			final Observation actual=new Observation(), expected=new Observation();
			reused.measure(q, r, self, actual);
			metric.newWorker().measure(q, r, self, expected);
			compare(expected, actual, 1);
		}
	}

	/** Different matrix workers repeatedly score the same inputs concurrently to expose shared selection. */
	private static void parallel(ProteinSequence query, ProteinSequence reference, FamilyAlignmentMetric custom,
		FamilyAlignmentMetric ordinary)throws Exception{
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread[] workers=new Thread[2];
		for(int i=0; i<2; i++){
			final boolean changed=i==0;
			final FamilyAlignmentMetric metric=changed ? custom : ordinary;
			workers[i]=new Thread(()->{
				try{
					final FamilyAlignmentMetric.Worker worker=metric.newWorker();
					final Observation out=new Observation();
					for(int j=0; j<100; j++){
						worker.measure(query, reference, 17, out);
						check(out.rawScore==(changed ? 64 : 12) && out.queryStart==(changed ? 0 : 1),
							"Concurrent matrix workers must not alter each other's optimum");
					}
				}catch(Throwable t){failure.compareAndSet(null, t);}
			});
			workers[i].start();
		}
		for(Thread worker:workers){worker.join();}
		if(failure.get()!=null){throw new AssertionError("Matrix worker isolation failed", failure.get());}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires identical states/coordinates and normalized metrics with exact integral unit scaling. */
	private static void compare(Observation x, Observation y, int scale){
		check(x.state==y.state && y.rawScore==scale*x.rawScore && y.denominatorSelf==scale*x.denominatorSelf &&
			y.fullQuerySelf==scale*x.fullQuerySelf && x.queryStart==y.queryStart && x.queryEnd==y.queryEnd &&
			x.referenceStart==y.referenceStart && x.referenceEnd==y.referenceEnd && x.columns==y.columns &&
				x.identities==y.identities,
			"Unit scaling must preserve state/geometry/identity and scale both score and denominator");
		if(x.measured()){
			check(x.identityPct()==y.identityPct() && x.blosum()==y.blosum(),
				"Unit scaling must preserve normalized metrics exactly");
		}
	}

	/** Literal published-score multiplication is the units-only control, not a new fitted matrix. */
	private static int[][] scaledScores(int scale){
		check(scale>0, "Fixture scale must be positive");
		final int[][] out=new int[20][20];
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){out[a][b]=scale*Blosum62.score((byte)a, (byte)b);}
		}
		return out;
	}

	/** Emits an explicit letter-ordered fixture so loader tests do not assume numeric encoding order. */
	private static String resource(ProteinSubstitutionMatrix matrix){
		check(matrix!=null, "Fixture matrix required");
		final StringBuilder out=new StringBuilder("#bbgap_substitution_matrix\t1\nname\t").append(matrix.name)
			.append("\nscale\t").append(matrix.scale).append("\ngap\t").append(matrix.linearGap).append('\n');
		for(char a:"ARNDCQEGHILKMFPSTWYV".toCharArray()){
			for(char b:"ARNDCQEGHILKMFPSTWYV".toCharArray()){
				out.append("score\t").append(a).append('\t').append(b).append('\t').append(matrix.score(
					(byte)AminoAcid.acidToNumberExtended[a], (byte)AminoAcid.acidToNumberExtended[b])).append('\n');
			}
		}
		return out.toString();
	}

	/** Creates one ASCII evidence file without overwriting and records its relative normalized path. */
	private static void write(Path path, String text)throws Exception{
		check(text!=null, "Fixture bytes required");
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		evidence.add(evidenceRoot.relativize(path.toAbsolutePath().normalize()).toString());
	}

	/** Creates a literal protein whose ID and ASCII sequence use the same supplied text. */
	private static ProteinSequence sequence(String s){
		return new ProteinSequence(s, s.getBytes(StandardCharsets.US_ASCII));
	}
	/** Enumerate created files directly; no directory discovery and no provenance hashing. */
	private static void finish(String marker)throws Exception{
		write(evidenceRoot.resolve("PASS"), marker+"\n");
		evidence.add("evidence.paths");
		evidence.add("evidence.links.tsv");
		Files.write(evidenceRoot.resolve("evidence.links.tsv"), new byte[0], StandardOpenOption.CREATE_NEW);
		Files.write(evidenceRoot.resolve("evidence.paths"), evidence, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		System.out.println(marker);
	}

	/** Absolute normalized root for the active fixture invocation. */
	private static Path evidenceRoot;
	/** Sorted unique generated-file paths, cleared when main begins. */
	private static final java.util.TreeSet<String> evidence=new java.util.TreeSet<String>();
	/** Requires IllegalArgumentException from the supplied malformed operation; other failures escape. */
	private static void expectFailure(Runnable action, String why){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException e){failed=true;}
		check(failed, why);
	}

	/** A rejection from an earlier unrelated guard is not evidence for the intended check. */
	private static void expectFailure(Runnable action, String diagnostic, String why){
		try{action.run();}catch(IllegalArgumentException e){
			check(e.getMessage()!=null && e.getMessage().contains(diagnostic), why+": wrong diagnostic: "+e.getMessage());
			return;
		}
		throw new AssertionError(why);
	}

	/** Throws AssertionError on a failed expectation, independent of -ea. */
	private static void check(boolean valid, String why){
		if(!valid){throw new AssertionError(why);}
	}
}
