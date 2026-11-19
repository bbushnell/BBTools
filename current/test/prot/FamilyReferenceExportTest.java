package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.LineParser1;
import prot.DigestSuffix;
import prot.FamilyCalibration;
import prot.FamilyCompactCodec;
import prot.FamilyCompactLoader;
import prot.FamilyCompactResourceValidate;
import prot.FamilyGateExport;
import prot.FamilyGateStats;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceLoader;
import prot.FamilyReferenceExport;
import prot.FamilyMemberScalars;
import prot.FamilyMemberPreparation;
import prot.FamilyTailCascade;
import prot.ProteinSequence;
import prot.ProteinSubstitutionMatrix;
import prot.ProteinSpectrum;
import shared.Shared;
import ddl.DDLFamilySearch;
import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import ddl.DDLPentamerBlacklist;
import structures.ByteBuilder;

/** Synthetic regression fixtures for family-reference export and compact resources.
 * Exercises conversion, directed rounding, scalar bounds, sparse/complete upper
 * evidence, monotone overrides, subsetting and resource binding. Search and
 * preparation cases invoke the real ranking/alignment code on small synthetic inputs;
 * this is not a biological accuracy or performance benchmark.
 * Run once in a dedicated JVM on Dori with fixture=NEW_DIR. The suite retains its
 * created files for inspection. Compression globals are restored around the explicit
 * BGZF checks; the suite does not promise restoration of every dependency global.
 * @author Sayu, Keqing, Fischl
 */
public final class FamilyReferenceExportTest{
	/*--------------------------------------------------------------*/
	/*----------------         Entry point          ----------------*/
	/*--------------------------------------------------------------*/

	/** Create the absent fixture directory and run all synthetic cases serially.
	 * Checks and expected rejections accumulate in the process-local counter; unexpected
	 * exceptions or failed checks abort immediately. The final suite success line is printed only at the end. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("fixture=")){
			throw new IllegalArgumentException("Require fixture=NEW_DIR");
		}
		final Path root=Paths.get(args[0].substring(8)).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final Fixture good=new Fixture(root.resolve("good"));
		final Path out=good.root.resolve("out");
		FamilyReferenceExport.run(good.args(out));
		check(read(out.resolve("PASS")).equals(FamilyReferenceExport.PASS+"\n"), "Export PASS");
		final ByteFile in=ByteFile.makeByteFile(out.resolve("references.tsv").toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[][] header=new byte[5][];
			for(int i=0; i<5; i++){header[i]=in.nextLine();}
			final FamilyReferenceCodec.Header h=FamilyReferenceCodec.parseHeader(header, 0, lp);
			check(h.totalFamilies==2 && h.liveFamilies==1 && h.emptyFamilies==1, "Full census");
			final byte[][] record=new byte[5][];
			for(int i=0; i<5; i++){record[i]=in.nextLine();}
			final FamilyReference f=FamilyReferenceCodec.parseFiveLines(record, 0, lp);
			check(f.id.equals("f0") && f.rank==0 && f.ddl[0]==1 && f.spectrum[0]==1 && f.spectrumNorm==1,
				"Reference data conversion");
			check(f.blosumFinal>=0.004 && Math.nextDown(f.blosumFinal)<0.004, "Original-double additive final gate");
			check(f.cosineAccept>=0.72 && Math.nextDown(f.cosineAccept)<0.72, "Original-double finite cosine floor");
			check(f.identityLower==30f && f.identityUpper==88f,
				"Unchanged member lower and once-applied1.1 identity margin survive packed export");
			for(int i=0; i<5; i++){record[i]=in.nextLine();}
			check(FamilyReferenceCodec.parseFiveLines(record, 0, lp).members==0 && in.nextLine()==null,
				"Empty family included, no extra records");
		}finally{check(!in.close(), "Reader success");}
		int cosinePair=0;
		for(double[] pair:new double[][]{{0.3, 0.3},{0.8, 0.2},{0.5, 0.96},{0.5, 0.9090909},{1, 0}}){
			final Fixture f=new Fixture(root.resolve("cosine_"+(cosinePair++)));
			cell(f.lower.resolve("cutoffs.tsv"), 1, 4, Double.toString(pair[0]));
			cell(f.upper.resolve("uppers.tsv"), 1, 3, Double.toString(pair[1]));
			f.bind();
			final Path target=f.root.resolve("out");
			FamilyReferenceExport.run(f.args(target));
			final List<String> rows=Files.readAllLines(target.resolve("references.tsv"), StandardCharsets.US_ASCII);
			final byte[][] record=new byte[5][];
			for(int i=0; i<5; i++){record[i]=rows.get(i+5).getBytes(StandardCharsets.US_ASCII);}
			final FamilyReference parsed=FamilyReferenceCodec.parseFiveLines(record, 0, new LineParser1('\t'));
			final double bound=Math.max(pair[0], Math.max(0.72, Math.max(pair[1]*1.1, pair[1]+0.05)));
			check(Float.isFinite(parsed.cosineAccept) && parsed.cosineAccept>=bound &&
				(double)Math.nextDown(parsed.cosineAccept)<bound,
				"Stored cosine must apply margin to original doubles once, including collapsed and >1 cases");
		}
		final Path again=good.root.resolve("again");
		FamilyReferenceExport.run(good.args(again));
		for(String name:new String[]{"references.tsv", "recipe.txt", "scalar_recipe.txt", "gate_quantiles.tsv",
			"inputs.tsv", "PASS", "result.sha80.tsv"}){
			check(Arrays.equals(Files.readAllBytes(out.resolve(name)), Files.readAllBytes(again.resolve(name))),
				"Deterministic export: "+name);
		}
		testCompactBinding(root);
		testSupplemental(good, out);
		testAlignmentReplacement(root);
		testMillionQueryBootstrap(root);
		testConsensusLengthContract(root);
		testHistoricalFloors(good, out);
		testScalars(good, out);
		testScalarApi(good, out);
		testPreparedReferences(root);
		testRaise(good, out);
		testMonomer(good.root, out);
		testSubset(good.root, out);
		testAdditional(root);
		reject(good.args(out));
		check(Files.exists(out.resolve("PASS")), "Existing output preserved");
		for(String mutation:new String[]{"seed", "profile_overflow", "duplicate_profile", "wrong_profile_count",
			"wrong_spectrum", "missing_upper",
			"duplicate_rank", "empty_metadata", "wrong_summary", "wrong_consensus", "duplicate_consensus", "bad_pin",
			"row_limit", "legacy_gate"}){
			final Fixture f=new Fixture(root.resolve(mutation));
			switch(mutation){
				case "seed":
					cell(f.refs.resolve("summary.tsv"), 1, 2, "1");
					break;
				case "profile_overflow":
					cell(f.refs.resolve("profiles.tsv"), 2, 7, array(65536));
					break;
				case "duplicate_profile":
					write(f.refs.resolve("profiles.tsv"),
						read(f.refs.resolve("profiles.tsv"))+line(f.refs.resolve("profiles.tsv"), 2)+"\n");
					break;
				case "wrong_profile_count":
					cell(f.refs.resolve("profiles.tsv"), 2, 3, "2");
					break;
				case "wrong_spectrum":
					cell(f.refs.resolve("spectra.tsv"), 1, 1, "unknown");
					break;
				case "missing_upper":
					write(f.upper.resolve("uppers.tsv"), FamilyGateExport.RAW_UPPER_HEADER+"\n");
					break;
				case "duplicate_rank":
					cell(f.refs.resolve("members.tsv"), 2, 0, "0");
					break;
				case "empty_metadata":
					cell(f.refs.resolve("members.tsv"), 2, 4, "1");
					break;
				case "wrong_summary":
					cell(f.refs.resolve("summary.tsv"), 1, 8, "2");
					break;
				case "wrong_consensus":
					write(f.consensus, ">unknown\nAAAA\n>empty\nX\n");
					break;
				case "duplicate_consensus":
					write(f.consensus, ">f0\nAAAA\n>f0\nX\n");
					break;
				case "legacy_gate":
					cell(f.upper.resolve("uppers.tsv"), 0, 6, "blosum_a_sq_final_a48f32_v2");
					break;
				default:
					break;
			}
			f.bind();
			final Path badOut=f.root.resolve("out");
			final String[] bad=f.args(badOut);
			if(mutation.equals("bad_pin")){replace(bad, "consensussha80", "00000000000000000000");}
			if(mutation.equals("row_limit")){replace(bad, "maxfamilies", "1");}
			reject(bad);
			check(!Files.exists(badOut), "Invalid input must not publish output: "+mutation);
		}
		final ProteinSpectrum.Spectrum s=new ProteinSpectrum.Spectrum(2);
		check(FamilyReferenceExport.frequencies(s)[0]==0, "Zero-total dimer normalization");
		s.counts[0]=Long.MAX_VALUE-1;
		s.counts[1]=1;
		s.total=Long.MAX_VALUE;
		check(FamilyReferenceExport.frequencies(s)[1]==(float)(1.0/Long.MAX_VALUE), "Divide in double before float cast");
		System.out.println("FAMILY_REFERENCE_EXPORT_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------  Projection and subsetting   ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal row sums distinguish sequential float addition from double accumulation and renormalization. */
	private static void testMonomer(Path root, Path baseline) throws Exception{
		final FamilyReference f=firstFamily(baseline);
		final float[] spectrum=new float[400];
		spectrum[0]=0.5f;
		for(int i=1; i<=16; i++){spectrum[i]=Math.scalb(1f, -26);}
		spectrum[20]=0.5f-Math.scalb(1f, -22);
		// Each tiny addend rounds away after0.5f; a double sum or normalized projection would differ.
		final float[] expected=new float[20];
		expected[0]=0.5f;
		expected[1]=spectrum[20];
		final double[] bounds={-100, 100, -100, 100, -100, 100, -100, 100, -100, 100};
		final char[] sketch=new char[1024];
		DDLFamilyProfiles.sketch(f.consensusRaw, 5, 0, sketch);
		final FamilyReference scalar=new FamilyReference(f.id, 0, 30, f.residues, f.validWords, f.consensus, sketch, spectrum,
			f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept, f.identityUpper, f.blosumFinal,
			f.consensusRaw, 31, bounds);
		final FamilyReference zero=new FamilyReference("zero", 1, 1, 100, 0, new ProteinSequence("zero", f.consensusRaw),
			new char[1024], new float[400], f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept,
			f.identityUpper, f.blosumFinal, f.consensusRaw);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA, 5, 1024, 0, 3, 2, 1,
			FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906);
		final ByteBuilder parentBytes=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		FamilyReferenceCodec.appendFiveLines(scalar, parentBytes);
		FamilyReferenceCodec.appendFiveLines(zero, parentBytes);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty", 2), parentBytes);
		final Path parent=root.resolve("monomer_parent.tsv");
		write(parent, parentBytes.toString());
		Path first=null;
		for(int threads:new int[]{1, 4}){
			final Path out=root.resolve("monomer_"+threads);
			FamilyReferenceExport.run(monomerArgs(parent, out, threads));
			check(read(out.resolve("PASS")).equals(FamilyReferenceExport.MONOMER_PASS+"\n"),
				"Compact conversion publishes its own PASS");
			final Path reference=out.resolve("references.tsv");
			final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(reference.toString(), threads, 3);
			check(loaded.header.spectrumWidth==20 && loaded.header.schema.equals(FamilyReferenceCodec.MONOMER_SCHEMA) &&
				loaded.header.scalarZ==3.8906, "Conversion declares20-value storage and retains scalar binding");
			for(int i=0; i<20; i++){
				check(Float.floatToRawIntBits(loaded.all[0].spectrum[i])==Float.floatToRawIntBits(expected[i]),
					"Projection retains exact sequential float row sum at "+i);
				check(loaded.all[1].spectrum[i]==0, "Live zero spectrum remains zero");
			}
			check(loaded.all[2].spectrum==null && loaded.live.length==2,
				"Empty roster entries remain distinct from live zero spectra");
			final List<String> before=Files.readAllLines(parent, StandardCharsets.US_ASCII),
				after=Files.readAllLines(reference, StandardCharsets.US_ASCII);
			check(before.size()==20 && after.size()==before.size(), "Conversion retains header plus three five-line records");
			for(int line=0; line<before.size(); line++){
				if(line!=0 && line!=4 && line!=9 && line!=14){
					check(before.get(line).equals(after.get(line)), "Non-spectrum bytes retained at line "+line);
				}
			}
			check(read(out.resolve("summary.tsv")).endsWith("\n3\t2\t1\t800\t40\n"), "Storage counts exclude empty families");
			check(read(out.resolve("inputs.tsv")).contains(parent+"\t"+hash(parent)), "Conversion records its pinned parent");
			for(String row:Files.readAllLines(out.resolve("result.sha80.tsv"), StandardCharsets.US_ASCII).subList(1, 6)){
				final String[] fields=row.split("\t");
				check(fields.length==2 && hash(out.resolve(fields[0])).equals(fields[1]), "Converted artifact manifest pin");
			}
			if(first==null){
				first=reference;
			}else{
				check(Arrays.equals(Files.readAllBytes(first), Files.readAllBytes(reference)),
					"Conversion is independent of loader thread setting");
			}
		}
		check(Arrays.equals(parentBytes.toBytes(), Files.readAllBytes(parent)), "Conversion never changes its parent");
		final Path repeat=root.resolve("monomer_repeated");
		reject(monomerArgs(first, repeat, 1));
		check(!Files.exists(repeat), "Already-projected input fails before creating output");
		final String[] badPin=monomerArgs(parent, root.resolve("monomer_bad_pin"), 1);
		replace(badPin, "parentsha80", "00000000000000000000");
		reject(badPin);
		final String[] base=monomerArgs(parent, root.resolve("monomer_bad_option"), 1),
			badOption=Arrays.copyOf(base, base.length+1);
		badOption[base.length]="identityfloor=40";
		reject(badOption);
		check(!Files.exists(root.resolve("monomer_bad_pin")) && !Files.exists(root.resolve("monomer_bad_option")),
			"Invalid conversion publishes no output");
		reject(monomerArgs(parent, first.getParent(), 1));
		final Path overrides=root.resolve("monomer_empty_overrides.tsv"), raised=root.resolve("monomer_raise_empty");
		write(overrides, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n");
		FamilyReferenceExport.run(raiseArgs(first, overrides, raised));
		check(Arrays.equals(Files.readAllBytes(first), Files.readAllBytes(raised.resolve("references.tsv"))),
			"Empty raise preserves compact schema and every byte");
		testCompactSearch(root, parent, first, f.consensusRaw);
	}
	/** Build fresh pinned monomer-conversion arguments with a ten-family limit. */
	private static String[] monomerArgs(Path parent, Path out, int threads){
		return new String[]{"mode=monomer", "parent="+parent, "parentsha80="+hash(parent), "out="+out, "maxfamilies=10",
			"threads="+threads};
	}

	/** Survivor subsetting removes retired/empty families before positional rows can enter the shortlist index. */
	private static void testSubset(Path root, Path baseline) throws Exception{
		final FamilyReference f=firstFamily(baseline);
		final FamilyReference[] records=new FamilyReference[4];
		for(int i=0; i<records.length; i++){
			records[i]=i==2 ? FamilyReference.empty("empty", i) :
				new FamilyReference("f"+i, i, 30, f.residues, f.validWords, new ProteinSequence("f"+i, f.consensusRaw),
					f.ddl, f.spectrum, f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept,
					f.identityUpper, f.blosumFinal, f.consensusRaw, f.scalarMask, f.scalarBounds);
		}
		final Path parent=root.resolve("subset_parent.tsv");
		writeFloor(parent, false, records);
		final Path source=root.resolve("subset_profiles.tsv"), summary=root.resolve("summary.tsv"),
			list=root.resolve("subset_list.tsv"),
			survivors=root.resolve("survivors.tsv"), blacklist=root.resolve("subset_blacklist.faa"),
			out=root.resolve("subset_out");
		write(blacklist, ">word\nCCCCC\n");
		final DDLPentamerBlacklist excluded=DDLPentamerBlacklist.load(blacklist, hash(blacklist));
		final ByteBuilder profiles=new ByteBuilder(DDLFamilyProfiles.HEADER16).nl();
		for(String id:new String[]{"f0", "f1", "f3"}){
			final String row=id.equals("f3") ? "\t1024\t31\t111\t7\t" : "\t1024\t30\t100\t1\t";
			profiles.append("tail\t").append(id).append(row).append("literal\t").append(array(1)).tab().append(ones())
				.tab().append(array(1)).nl();
			profiles.append("tail\t").append(id).append(row).append("nonempty\t").append(array(1)).tab()
				.append(array(1)).tab().append(array(1)).nl();
		}
		write(source, profiles.toString());
		write(summary, "namespace\tmanifest_sha80\tfamilies_sha80\tchunks_sha80\tk\tseed\texponent\tregisterbits\tbuckets\tshard\tshards"+
			"\tmanifest_families\tbuilt_families\tzero_families\tmembers\tresidues\tvalid_words\tload_seconds\tbuild_seconds\tsketch"+
			"\tpartitions\toverlap\ttag_bits\tsegment_words\tblacklist_sha80\tblacklist_words\n"+
			"tail\t00000000000000000000\t11111111111111111111\t22222222222222222222\t5\t0\t6\t16\t1024\t0\t1"+
			"\t3\t3\t0\t91\t311\t9\t3\t7\tpositional-v1\t4\t2\t2\t49\t"+excluded.sha80()+"\t"+excluded.size()+"\n");
		write(list, source+"\t"+hash(source)+"\t"+hash(summary)+"\n");
		write(survivors, "rank\tfamily_id\n0\tf0\n3\tf3\n");
		FamilyReferenceExport.run(subsetArgs(parent, survivors, list, out));
		check(read(out.resolve("PASS")).equals(FamilyReferenceExport.SUBSET_PASS+"\n"), "Subset publishes its own PASS");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(out.resolve("references.tsv").toString(), 2, 10);
		check(loaded.header.totalFamilies==2 && loaded.header.liveFamilies==2 && loaded.header.emptyFamilies==0 &&
			loaded.all[0].id.equals("f0") && loaded.all[0].rank==0 && loaded.all[1].id.equals("f3") &&
			loaded.all[1].rank==1, "Subset densifies only surviving live rows in parent order");
		final List<String> parentLines=Files.readAllLines(parent, StandardCharsets.US_ASCII),
			subsetLines=Files.readAllLines(out.resolve("references.tsv"), StandardCharsets.US_ASCII);
		for(int i=1; i<5; i++){
			check(parentLines.get(5+i).equals(subsetLines.get(5+i)) && parentLines.get(20+i).equals(subsetLines.get(10+i)),
				"Subset preserves survivor payload bytes at line "+i);
		}
		final DDLFamilySearch.Index index=DDLFamilySearch.loadForAttribution(out.resolve("ddlprofiles.tsv"), 1024, 5, 0,
			"positional-v1", 4, 2, excluded);
		check(index.families()==2 && index.names[0].equals("tail:f0") && index.names[1].equals("tail:f3"),
			"Subset profile list loads with retired families physically absent");
		check(read(out.resolve("survivor_roster.tsv")).contains("3\tf3\t1\t30\n") &&
			read(out.resolve("summary.tsv")).contains("\t2\t2\t0\t61\t211\t8\t0\t0\tpositional-v1\t4\t2\t2\t49\t"+excluded.sha80()+"\t1\n"),
			"Subset preserves reference counts while rewriting summary counts from profile rows");
		for(String mutation:new String[]{"duplicate", "missing_profile", "empty"}){
			final Path badRoot=root.resolve("subset_bad_"+mutation), badSurvivors=badRoot.resolve("survivors.tsv"),
				badList=badRoot.resolve("list.tsv");
			Files.createDirectory(badRoot);
			write(badSurvivors, mutation.equals("duplicate") ? "family_id\nf0\nf0\n" :
				mutation.equals("empty") ? "family_id\nempty\n" : "family_id\nf0\nf3\n");
			if(mutation.equals("missing_profile")){
				final Path missing=badRoot.resolve("profiles.tsv"), badSummary=badRoot.resolve("summary.tsv");
				final String row="\t1024\t30\t100\t1\t";
				write(missing, DDLFamilyProfiles.HEADER16+"\n"+
					"tail\tf0"+row+"literal\t"+array(1)+"\t"+ones()+"\t"+array(1)+"\n"+
					"tail\tf0"+row+"nonempty\t"+array(1)+"\t"+array(1)+"\t"+array(1)+"\n");
				write(badSummary, read(summary));
				write(badList, missing+"\t"+hash(missing)+"\t"+hash(badSummary)+"\n");
			}
			reject(subsetArgs(parent, badSurvivors, mutation.equals("missing_profile") ? badList : list, badRoot.resolve("out")));
		}
	}

	/** Build fresh pinned survivor/profile arguments with the fixed synthetic loader settings. */
	private static String[] subsetArgs(Path parent, Path survivors, Path profiles, Path out){
		return new String[]{"mode=subset", "parent="+parent, "parentsha80="+hash(parent), "survivors="+survivors,
			"survivorssha80="+hash(survivors), "ddlprofiles="+profiles, "ddlprofilessha80="+hash(profiles),
			"out="+out, "maxfamilies=10", "threads=2"};
	}

	/*--------------------------------------------------------------*/
	/*----------------  Compact resource fixtures   ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercise compact rendering, direct BGZF output, resource binding and strict sidecar checks.
	 * Uses one live and one empty synthetic family, literal header expectations,
	 * byte equality and deliberately mismatched resources. No production resource is changed. */
	private static void testCompactBinding(Path root) throws Exception{
		final Path in=Files.createDirectory(root.resolve("compact_inputs"));
		final String sourceProvenance="original_component_v1";
		compactInputs(in, sourceProvenance);
		final Path alphabet=in.resolve("alphabet.config"), bounds=in.resolve("scalar_bounds.tsv"),
			blacklist=in.resolve("source_blacklist.faa"), neutralBlacklist=in.resolve("blacklist.faa"),
			matrix=in.resolve("possum55.tsv"), out=root.resolve("compact_out");
		final String aabits=dimerAabits(120);
		check(dimerGroups(aabits)==120, "Synthetic compact alphabet fixture has 120 groups");
		write(alphabet, "#bbgap_alphabet_config\t1\naabits="+aabits+"\n");
		blacklist(blacklist, aabits, sourceProvenance);
		blacklist(neutralBlacklist, aabits);
		final String binding=possum55(matrix);
		write(bounds, "scalar\trows\tmin\tmax\tmin_family\tmax_family\n"+
			"hydrophobicity\t79565\t-2.740000\t1.777679\ttail_75534\ttail_48259\n"+
			"net_charge\t79565\t-0.621212\t0.481675\ttail_48392\ttail_72959\n"+
			"aromaticity\t79565\t0\t0.290000\ttail_5426\ttail_62434\n"+
			"fraction_small\t79565\t0.030000\t0.917526\ttail_69058\ttail_38096\n"+
			"pI_tail\t79565\t1.9168\t13.4343\ttail_70514\ttail_64756\n");
		runCompactFixture(compactArgs(in, alphabet, bounds, blacklist, neutralBlacklist, matrix, binding, out));
		check(read(out.resolve("PASS")).equals(FamilyReferenceExport.COMPACT_PASS+"\n"),
			"Compact binding publishes its own PASS");
		final List<String> rows=Files.readAllLines(out.resolve("families.tsv"), StandardCharsets.US_ASCII);
		final byte[][] header=new byte[10][];
		for(int i=0; i<10; i++){header[i]=rows.get(i).getBytes(StandardCharsets.US_ASCII);}
		final FamilyCompactCodec.Header h=FamilyCompactCodec.parseHeader(header, 0, new LineParser1('\t'));
		check(h.get("view").equals("fixture") && h.get("families").equals("2") && h.get("live").equals("1") &&
			h.get("empty").equals("1"), "Compact header binds the synthetic census");
		check(h.get("alphabet").equals("alphabet.config") && h.get("blacklist").equals("blacklist.faa.gz") &&
			h.get("matrixresource").equals("POSSUM55.tsv") && h.get("matrix").equals("POSSUM55"),
			"Compact header stores package-root shared resource locators");
		check(h.get("geometry").startsWith("dimer6+1,3,5+") && h.get("matrixdisplay").equals("POSSUM55") &&
			h.get("packageprovenance").equals("bbgap_release_20261116_v1") &&
			h.get("blacklistprovenance").equals("bbgap_blacklist_k6_2048_20261116_v1") &&
			h.get("hydrophobicity_min").equals("-2.740000") && h.get("aromaticity_max").equals("0.290000") &&
			h.get("fraction_small_min").equals("0.030000") && h.get("pI_tail_max").equals("13.4343"),
			"Compact header binds accepted Q76 geometry and scalar resource metadata");
		FamilyCompactLoader.load(out.resolve("families.tsv").toString(), 2, 2);
		check(read(out.resolve("recipe.txt")).contains("payload_rendering=complete\n") &&
			read(out.resolve("recipe.txt")).contains("neutralblacklist="+neutralBlacklist+"\n") &&
			read(out.resolve("recipe.txt")).contains("output=families.tsv\n") &&
			Files.readAllLines(out.resolve("inputs.tsv"), StandardCharsets.US_ASCII).size()==14,
			"Compact binding records structural inputs and renders families.tsv");
		final Path directOut=root.resolve("compact_direct"), directDir=root.resolve("compact_direct_payload"),
			directFamilies=directDir.resolve("families.tsv.gz");
		Files.createDirectory(directDir);
		final String[] base=compactArgs(in, alphabet, bounds, blacklist, neutralBlacklist, matrix, binding, directOut),
			direct=Arrays.copyOf(base, base.length+2);
		direct[base.length]="outfamilies="+directFamilies;
		direct[base.length+1]="threads=2";
		final TestZipState saved=TestZipState.capture();
		try{
			TestZipState.setNondefault();
			final TestZipState expected=TestZipState.capture();
			runCompactFixture(direct);
			check(expected.equalsCurrent(), "Explicit compact BGZF restores zip globals after success");
			exerciseCompactBgzipFailure(root.resolve("compact_direct_failure.tsv.gz"));
			check(expected.equalsCurrent(), "Explicit compact BGZF restores zip globals after writer failure");
		}finally{saved.restore();}
		check(!Files.exists(directOut.resolve("families.tsv")) &&
			readCompressed(directFamilies).equals(read(out.resolve("families.tsv"))),
			"Explicit compact BGZF output preserves uncompressed payload bytes without a raw family file");
		FamilyCompactLoader.load(directFamilies.toString(), 2, 2);
		exerciseCompactDecimalSidecars(root, in, alphabet, bounds, blacklist, neutralBlacklist, matrix, binding);
		final Path acceptedMatrix=in.resolve("custom40.tsv"), badAcceptedMatrix=in.resolve("custom40_bad.tsv");
		matrix(acceptedMatrix, "custom40", 20, -1);
		matrix(badAcceptedMatrix, "custom40", 19, -1);
		FamilyCompactResourceValidate.main(new String[]{"mode=matrixbinding", "matrixresource="+matrix,
			"acceptedmatrix="+acceptedMatrix});
		reject(()->FamilyCompactResourceValidate.main(new String[]{"mode=matrixbinding",
			"matrixresource="+matrix, "acceptedmatrix="+badAcceptedMatrix}));
		FamilyCompactResourceValidate.main(new String[]{"mode=blacklist", "short="+blacklist,
			"tails="+blacklist, "out="+root.resolve("blacklist_neutral.faa")});
		final Path badPolicyBlacklist=in.resolve("blacklist_bad.faa");
		blacklist(badPolicyBlacklist, dimerAabits(119));
		reject(()->FamilyCompactResourceValidate.main(new String[]{"mode=blacklist", "short="+blacklist,
			"tails="+badPolicyBlacklist, "out="+root.resolve("blacklist_bad_out.faa")}));
		final Path badBlacklistKey=in.resolve("blacklist_bad_key.faa");
		blacklist(badBlacklistKey, aabits, 0, absentBlacklistWord(aabits));
		final Path badKeyOut=root.resolve("blacklist_bad_key_out.faa");
		rejectMessage(()->FamilyCompactResourceValidate.main(new String[]{"mode=blacklist", "short="+blacklist,
			"tails="+badBlacklistKey, "out="+badKeyOut}), "loaded keys differ");
		check(!Files.exists(badKeyOut), "Loaded-key blacklist rejection must not write a neutral output");
		reject(()->runCompactFixture(compactArgs(in, alphabet, bounds, neutralBlacklist, neutralBlacklist, matrix, binding,
			root.resolve("compact_bad_source_provenance"))));
		reject(()->runCompactFixture(compactArgs(in, alphabet, bounds, blacklist, badBlacklistKey, matrix, binding,
			root.resolve("compact_bad_neutral_keys"))));
		final Path validate=root.resolve("compact_validate");
		runCompactValidate(validate, directFamilies, in, alphabet, bounds, binding);
		check(read(validate.resolve("PASS")).contains("view=fixture physical=2 live=1 members=30"),
			"Compact validator accepts the synthetic direct-BGZF payload");
		final String badBinding=binding.replace(":20:", ":19:");
		check(!badBinding.equals(binding), "Matrix-header negative must mutate a score field");
		reject(()->runCompactValidate(root.resolve("compact_validate_bad_header"), directFamilies, in, alphabet, bounds,
			badBinding));
		final Path badHydro=in.resolve("hydro_bad.tsv");
		write(badHydro, "family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t31\t0.5\n1\tempty\t0\tNA\n");
		reject(()->runCompactValidate(root.resolve("compact_validate_bad_hydro"), directFamilies, in, alphabet, bounds,
			binding, badHydro, in.resolve("profiles.tsv")));
		final Path badProfiles=in.resolve("profiles_bad.tsv");
		final String profiles=read(in.resolve("profiles.tsv"));
		final String mutatedProfiles=profiles.replace("\tnonempty\t1,1,1", "\tnonempty\t2,1,1");
		check(!mutatedProfiles.equals(profiles), "Profile-vector negative must mutate the live DDL oracle");
		write(badProfiles, mutatedProfiles);
		reject(()->runCompactValidate(root.resolve("compact_validate_bad_profiles"), directFamilies, in, alphabet, bounds,
			binding, in.resolve("hydro.tsv"), badProfiles));
		final String[] bad=compactArgs(in, alphabet, bounds, blacklist, neutralBlacklist, matrix, binding,
			root.resolve("compact_bad_locator"));
		final String[] badLocator=Arrays.copyOf(bad, bad.length+1);
		badLocator[bad.length]="blacklistlocator=../blacklist.faa.gz";
		reject(()->runCompactFixture(badLocator));
		check(!Files.exists(root.resolve("compact_bad_locator")), "Escaping compact locators fail before output");
		for(int groups:new int[]{119, 121}){
			final Path badAlphabet=in.resolve("alphabet_"+groups+".config"),
				badBlacklist=in.resolve("blacklist_"+groups+".faa");
			final String badAabits=dimerAabits(groups);
			check(dimerGroups(badAabits)==groups, "Synthetic compact negative alphabet has requested groups");
			write(badAlphabet, "#bbgap_alphabet_config\t1\naabits="+badAabits+"\n");
			blacklist(badBlacklist, badAabits);
			reject(()->runCompactFixture(compactArgs(in, badAlphabet, bounds, badBlacklist, badBlacklist, matrix, binding,
				root.resolve("compact_bad_groups_"+groups))));
			check(!Files.exists(root.resolve("compact_bad_groups_"+groups)),
				"Non-Q70 compact dimer group count fails before output: "+groups);
		}
	}

	/** Validate the synthetic compact payload against all default fixture sidecars. */
	private static void runCompactValidate(Path out, Path compact, Path in, Path alphabet, Path bounds,
			String binding) throws Exception{
		runCompactValidate(out, compact, in, alphabet, bounds, binding, in.resolve("hydro.tsv"), in.resolve("profiles.tsv"));
	}

	/** Validate with substituted hydrophobicity and profile sidecars for rejection cases. */
	private static void runCompactValidate(Path out, Path compact, Path in, Path alphabet, Path bounds,
			String binding, Path hydro, Path profiles) throws Exception{
		runCompactValidate(out, compact, in, alphabet, bounds, binding, in.resolve("lengths.tsv"), hydro,
			in.resolve("family_gc.tsv"), in.resolve("family_map.tsv"), profiles);
	}

	/** Invoke the public resource validator with the two-family synthetic census and supplied sidecars. */
	private static void runCompactValidate(Path out, Path compact, Path in, Path alphabet, Path bounds,
			String binding, Path lengths, Path hydro, Path gc, Path gcmap, Path profiles) throws Exception{
		FamilyCompactResourceValidate.main(new String[]{"view=fixture", "references="+in.resolve("references.tsv"),
			"compact="+compact, "profilevectors="+profiles, "lengthmedians="+lengths,
			"familyhydro="+hydro, "familygc="+gc,
			"familygcmap="+gcmap, "scalarbounds="+bounds,
			"alphabet="+alphabet, "matrixbinding="+binding, "out="+out,
			"families=2", "live=1", "empty=1", "members=30", "minaa=40", "maxaa=99",
			"alignmentminspan=28", "alignmentminpaired=20", "marginscale=0.985", "threads=2"});
	}

	/** Invoke the private compact exporter with an explicit small census via reflection.
	 * This avoids changing production census rules. Unwrap the target exception/error
	 * so negative cases observe the exporter failure rather than the reflection wrapper. */
	private static void runCompactFixture(String[] args) throws Exception{
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String arg:args){
			final String[] pair=arg.split("=", 2);
			options.put(pair[0], pair[1]);
		}
		final Class<?> census=Class.forName("prot.FamilyReferenceExport$CompactCensus");
		final java.lang.reflect.Constructor<?> ctor=census.getDeclaredConstructor(String.class, int.class,
			int.class, int.class, long.class, int.class, int.class, int.class, int.class, String.class);
		ctor.setAccessible(true);
		final Object fixture=ctor.newInstance("fixture", 2, 1, 1, 30L, 40, 99, 28, 20, "0.985");
		final java.lang.reflect.Method compact=FamilyReferenceExport.class.getDeclaredMethod("compact",
			HashMap.class, census);
		compact.setAccessible(true);
		try{
			compact.invoke(null, options, fixture);
		}catch(java.lang.reflect.InvocationTargetException failure){
			final Throwable cause=failure.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Pass a null family to the private direct-BGZF writer and require NullPointerException.
	 * The enclosing caller checks global restoration; partial output is retained for inspection. */
	private static void exerciseCompactBgzipFailure(Path out) throws Exception{
		final java.lang.reflect.Method writer=FamilyReferenceExport.class.getDeclaredMethod("writeCompactFamilies",
			Path.class, boolean.class, int.class, ByteBuilder.class, FamilyReference[].class, int[].class,
			double[].class, double[].class, long[].class, double[].class, double[].class, char[][].class);
		writer.setAccessible(true);
		try{
			writer.invoke(null, new Object[]{out, true, 2, new ByteBuilder("partial\n"),
				new FamilyReference[]{null}, new int[1], new double[1], new double[1], new long[1],
				new double[1], new double[1], new char[1][0]});
		}catch(java.lang.reflect.InvocationTargetException failure){
			final Throwable cause=failure.getCause();
			if(cause instanceof NullPointerException){
				checks++;
				return;
			}
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
		throw new AssertionError("Malformed direct BGZF writer fixture accepted");
	}

	/** Check decimal/exponent sidecars against literal expected raw binary64 values.
	 * Render and independently reload/validate each case using the same synthetic payload. */
	private static void exerciseCompactDecimalSidecars(Path root, Path in, Path alphabet, Path bounds,
			Path blacklist, Path neutralBlacklist, Path matrix, String binding) throws Exception{
		final String[][] variants={
			{"length_exponent", "5.05e1", "0.5", "0.5", "0.1",
				"4049400000000000", "3fe0000000000000", "3fe0000000000000", "3fb999999999999a"},
			{"hydro_decimal", "50", "-0.3", "0.5", "0.1",
				"4049000000000000", "bfd3333333333333", "3fe0000000000000", "3fb999999999999a"},
			{"hydro_exponent", "50", "-3e-1", "0.5", "0.1",
				"4049000000000000", "bfd3333333333333", "3fe0000000000000", "3fb999999999999a"},
			{"gc_decimal", "50", "0.5", "0.3", "0.3",
				"4049000000000000", "3fe0000000000000", "3fd3333333333333", "3fd3333333333333"},
			{"gc_exponent", "50", "0.5", "3e-1", "3e-1",
				"4049000000000000", "3fe0000000000000", "3fd3333333333333", "3fd3333333333333"}};
		for(String[] v:variants){
			final Path sidecars=Files.createDirectory(root.resolve("compact_decimal_"+v[0])),
				lengths=sidecars.resolve("lengths.tsv"), hydro=sidecars.resolve("hydro.tsv"),
				gc=sidecars.resolve("family_gc.tsv"), map=sidecars.resolve("family_map.tsv"),
				out=root.resolve("compact_decimal_out_"+v[0]);
			write(lengths, "family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
				"f0\t30\t50\t"+v[1]+"\t51\t30\nempty\t0\tNA\tNA\tNA\t0\n");
			write(hydro, "family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t30\t"+v[2]+"\n1\tempty\t0\tNA\n");
			write(gc, "stable_family_id\troster_rank\tfamily_id\tmembers\tresidues\tmeasured\tmean_gc\tsample_sd\ttail_min_aa\ttail_max_aa\n"+
				"7\t0\tf0\t30\t1500\t30\t"+v[3]+"\t"+v[4]+"\t40\t99\n");
			write(map, "stable_family_id\troster_rank\tfamily_id\tassigned\tresidues\n7\t0\tf0\t30\t1500\n");
			final String[] args=compactArgs(in, alphabet, bounds, blacklist, neutralBlacklist, matrix, binding, out);
			replace(args, "lengthmedians", lengths.toString());
			replace(args, "familyhydro", hydro.toString());
			replace(args, "familygc", gc.toString());
			replace(args, "familygcmap", map.toString());
			runCompactFixture(args);
			final FamilyCompactLoader.Loaded loaded=FamilyCompactLoader.load(out.resolve("families.tsv").toString(), 2, 2);
			check(Double.doubleToRawLongBits(loaded.lengthMedians[0])==Long.parseUnsignedLong(v[5], 16) &&
				Double.doubleToRawLongBits(loaded.kdMeans[0])==Long.parseUnsignedLong(v[6], 16) &&
				Double.doubleToRawLongBits(loaded.gcMeans[0])==Long.parseUnsignedLong(v[7], 16) &&
				Double.doubleToRawLongBits(loaded.gcSampleSds[0])==Long.parseUnsignedLong(v[8], 16),
				"Compact sidecar decimals preserve exporter raw64 parsing: "+v[0]);
			runCompactValidate(root.resolve("compact_decimal_validate_"+v[0]), out.resolve("families.tsv"),
				in, alphabet, bounds, binding, lengths, hydro, gc, map, in.resolve("profiles.tsv"));
		}
	}

	/** Write one live/one empty legacy reference and matching compact-renderer sidecars.
	 * The 30-member live family uses 512 registers and a 50-residue consensus; all data are synthetic. */
	private static void compactInputs(Path in, String blacklistProvenance) throws Exception{
		final char[] ddl=new char[512];
		Arrays.fill(ddl, (char)1);
		final float[] monomers=new float[20];
		monomers[0]=1;
		final double[] scalars={-1, 1, -1, 1, 0, .1, .1, .9, 2, 10};
		final byte[] consensus="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes(StandardCharsets.US_ASCII);
		final FamilyReference live=new FamilyReference("f0", 0, 30, 1500, 150,
			new ProteinSequence("f0", consensus), ddl, monomers, .5f, 30, .1f, .75f, 40, .2f,
			consensus, 31, scalars);
		final FamilyReferenceCodec.Header h=new FamilyReferenceCodec.Header(FamilyReferenceCodec.MONOMER_SCHEMA,
			6, 512, 0, 2, 1, 1, FamilyReferenceCodec.FORMULA_VERSION,
			FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906);
		final ByteBuilder refs=FamilyReferenceCodec.appendHeader(h, new ByteBuilder());
		FamilyReferenceCodec.appendFiveLines(live, refs, h);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty", 1), refs, h);
		write(in.resolve("references.tsv"), refs.toString());
		write(in.resolve("ddlprofiles.tsv"), in.resolve("profiles.tsv")+"\tunmeasured\tunmeasured\n");
		write(in.resolve("profiles.tsv"), DDLFamilyProfiles.HEADER16+"\ntail\tf0\t512\t30\t1500\t150\tnonempty\t"+
			array(1, 512)+"\t"+array(30, 512)+"\t"+array(30, 512)+"\n");
		write(in.resolve("summary.tsv"), "namespace\tmanifest_sha80\tk\tseed\texponent\tmanifest_families\tbuilt_families\tmembers\tresidues\tvalid_words\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\tgeometry\tblacklist_sha80\tblacklist_words\tblacklist_k\tblacklist_mode\tblacklist_provenance\n"+
			"tail\tunmeasured\t6\t0\t6\t1\t1\t30\t1500\t150\tpositional-v1\t4\t2\t2\t49\t"+
			ProteinKmerGeometry.make(6, dimerAabits(120)).name()+
			"\tunmeasured\t2048\t6\treduced_literal\t"+blacklistProvenance+"\n");
		write(in.resolve("lengths.tsv"), "family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
			"f0\t30\t50\t50\t50\t30\nempty\t0\tNA\tNA\tNA\t0\n");
		write(in.resolve("hydro.tsv"), "family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t30\t0.5\n1\tempty\t0\tNA\n");
		write(in.resolve("family_gc.tsv"),
			"stable_family_id\troster_rank\tfamily_id\tmembers\tresidues\tmeasured\tmean_gc\tsample_sd\ttail_min_aa\ttail_max_aa\n"+
			"7\t0\tf0\t30\t1500\t30\t0.5\t0.1\t40\t99\n");
		write(in.resolve("family_map.tsv"), "stable_family_id\troster_rank\tfamily_id\tassigned\tresidues\n7\t0\tf0\t30\t1500\n");
	}

	/** Return width repeated decimal values separated by commas; callers provide a positive width. */
	private static String array(int value, int width){
		final StringBuilder out=new StringBuilder();
		for(int i=0; i<width; i++){
			if(i>0){out.append(',');}
			out.append(value);
		}
		return out.toString();
	}

	/** Build compact-export arguments naming sidecars and shared resources.
	 * The fixture invocation substitutes its small census for the production short view. */
	private static String[] compactArgs(Path in, Path alphabet, Path bounds, Path blacklist, Path neutralBlacklist,
			Path matrix, String binding, Path out){
		return new String[]{"mode=compact", "view=short", "out="+out,
			"references="+in.resolve("references.tsv"), "ddlprofiles="+in.resolve("ddlprofiles.tsv"),
			"profilevectors="+in.resolve("profiles.tsv"), "profilesummary="+in.resolve("summary.tsv"),
			"lengthmedians="+in.resolve("lengths.tsv"), "familyhydro="+in.resolve("hydro.tsv"),
			"familygc="+in.resolve("family_gc.tsv"), "familygcmap="+in.resolve("family_map.tsv"),
			"scalarbounds="+bounds, "alphabet="+alphabet, "blacklist="+blacklist,
			"neutralblacklist="+neutralBlacklist,
			"matrixresource="+matrix, "matrixbinding="+binding};
	}
	/** Partition all 400 ordered amino-acid pairs into the requested synthetic groups.
	 * The first groups minus one pairs are singleton groups; all remaining pairs share the last. */
	private static String dimerAabits(int groups){
		check(groups>1 && groups<=400, "Dimer fixture group count");
		final String amino="ACDEFGHIKLMNPQRSTVWY";
		final StringBuilder out=new StringBuilder("dimer:1,3,5:");
		int pair=0;
		for(int i=0; i<amino.length(); i++){
			for(int j=0; j<amino.length(); j++){
				if(pair>0){out.append(pair<groups ? '/' : ',');}
				out.append(amino.charAt(i)).append(amino.charAt(j));
				pair++;
			}
		}
		return out.toString();
	}
	/** Count slash-delimited groups in this fixture's dimer:1,3,5: representation. */
	private static int dimerGroups(String aabits){
		return aabits.substring("dimer:1,3,5:".length()).split("/", -1).length;
	}
	/** Write the synthetic 2,048-word blacklist with the neutral release provenance token. */
	private static void blacklist(Path path, String aabits) throws Exception{
		blacklist(path, aabits, "bbgap_blacklist_k6_2048_20261116_v1");
	}
	/** Write the synthetic blacklist with explicit source or neutral provenance. */
	private static void blacklist(Path path, String aabits, String provenance) throws Exception{
		blacklist(path, aabits, -1, null, provenance);
	}
	/** Write a neutral blacklist replacing one selected word to exercise loaded-key rejection. */
	private static void blacklist(Path path, String aabits, int replace, String replacement) throws Exception{
		blacklist(path, aabits, replace, replacement, "bbgap_blacklist_k6_2048_20261116_v1");
	}
	/** Write 2,048 deterministic six-residue representatives from three dimer groups.
	 * A negative replacement index leaves all words unchanged; a selected index uses
	 * the supplied replacement. This helper intentionally permits malformed test inputs. */
	private static void blacklist(Path path, String aabits, int replace, String replacement, String provenance) throws Exception{
		final String[] groups=aabits.substring("dimer:1,3,5:".length()).split("/");
		for(int i=0; i<groups.length; i++){groups[i]=groups[i].split(",")[0];}
		final ByteBuilder out=new ByteBuilder(">bbgap_blacklist_v1;k=6;words=2048;mode=reduced_literal;provenance=");
		out.append(provenance).append(";aabits=");
		out.append(aabits).nl();
		for(int i=0; i<2048; i++){
			if(i>0){out.append(">b_").append(i).nl();}
			if(i==replace){
				out.append(replacement);
			}else{
				out.append(groups[i%groups.length]).append(groups[(i/groups.length)%groups.length])
					.append(groups[i/(groups.length*groups.length)]);
			}
			out.nl();
		}
		write(path, out.toString());
	}
	/** Return a representative outside the first 2,048 words for the tested 120-group alphabet. */
	private static String absentBlacklistWord(String aabits){
		final String[] groups=aabits.substring("dimer:1,3,5:".length()).split("/");
		return groups[0].split(",")[0]+groups[groups.length-1].split(",")[0]+groups[0].split(",")[0];
	}
	/** Write a synthetic matrix named POSSUM55 and return its actual loader binding.
	 * Scores are the fixture's 20/-1 values, not a claim about the scientific POSSUM55 table. */
	private static String possum55(Path path) throws Exception{
		matrix(path, ProteinSubstitutionMatrix.POSSUM55_NAME, 20, -1);
		return ProteinSubstitutionMatrix.load(ProteinSubstitutionMatrix.POSSUM55_NAME, path.toString()).binding();
	}
	/** Write a complete synthetic 20-symbol substitution matrix with fixed scale/gap metadata.
	 * Caller-supplied diagonal and off-diagonal scores define the binding test. */
	private static void matrix(Path path, String name, int diagonal, int offDiagonal) throws Exception{
		final String letters="ARNDCEQGHILKMFPSTWYV";
		final ByteBuilder out=new ByteBuilder("#bbgap_substitution_matrix\t1\n");
		out.append("name\t").append(name).nl();
		if(ProteinSubstitutionMatrix.POSSUM55_NAME.equals(name)){
			out.append("description\t").append(ProteinSubstitutionMatrix.POSSUM55_DESCRIPTION).nl();
		}
		out.append("scale\t4\ngap\t16\n");
		for(int i=0; i<letters.length(); i++){
			for(int j=0; j<letters.length(); j++){
				out.append("score").tab().append(letters.charAt(i)).tab().append(letters.charAt(j))
					.tab().append(i==j ? diagonal : offDiagonal).nl();
			}
		}
		write(path, out.toString());
	}

	/** Exercise real ranking/DP plus early rejection; incompatible modes must fail before any writer opens. */
	private static void testCompactSearch(Path root, Path dimer, Path monomer, byte[] sequence) throws Exception{
		final Path query=root.resolve("compact_queries.faa"), bounds=root.resolve("compact_bounds.tsv");
		final char[] miss=new char[sequence.length];
		Arrays.fill(miss, 'D');
		write(query, ">hit\n"+new String(sequence, StandardCharsets.US_ASCII)+
			"\n>miss\n"+new String(miss)+"\n>zero\nXXXXXXXXXX\n");
		write(bounds, "scalar\trows\tmin\tmax\tmin_family\tmax_family\n"+
			"hydrophobicity\t79565\t-2.740000\t1.777679\tlo\thi\n"+
			"net_charge\t79565\t-0.621212\t0.481675\tlo\thi\n"+
			"aromaticity\t79565\t0\t0.290000\tlo\thi\n"+
			"fraction_small\t79565\t0.030000\t0.917526\tlo\thi\n"+
			"pI_tail\t79565\t1.9168\t13.4343\tlo\thi\n");
		byte[] expected=null, expectedAudit=null;
		for(int threads:new int[]{1, 4}){
			for(int compact=0; compact<2; compact++){
				final Path out=root.resolve("compact_search_"+threads+"_"+compact+".tsv"), audit=Paths.get(out+".audit");
				final String[] args=compactSearchArgs(compact==0 ? dimer : monomer, query, out, audit, bounds, threads);
				FamilyTailCascade.run(args);
				final byte[] actual=Files.readAllBytes(out), actualAudit=Files.readAllBytes(audit);
				check(read(out).contains("hit\tassigned\tf0\t") && read(audit).contains("hit\tidentity_accept\t"),
					"Fixture must reach alignment and accept its exact consensus, not only compare two empty outputs");
				check(Files.readAllLines(out, StandardCharsets.US_ASCII).size()==4 && read(out).contains("zero\tno_hit\t"),
					"Search retains all three queries including a zero-spectrum rejection");
				if(expected==null){
					expected=actual;
					expectedAudit=actualAudit;
				}else{
					check(Arrays.equals(expected, actual) && Arrays.equals(expectedAudit, actualAudit),
						"400/20 storage and thread settings preserve assignment/audit bytes");
				}
			}
		}
		for(String mode:new String[]{"dimer", "monomer", "centered_monomer", "full", "off", "rankout"}){
			for(int compact=0; compact<2; compact++){
				final Path out=root.resolve("compact_mode_"+mode+"_"+compact+".tsv"), audit=Paths.get(out+".audit"),
					rank=Paths.get(out+".rank");
				final java.util.ArrayList<String> args=new java.util.ArrayList<String>(Arrays.asList(
					compactSearchArgs(compact==0 ? dimer : monomer, query, out, audit, bounds, 1)));
				if(mode.equals("rankout")){
					args.add("rankout="+rank);
				}else{
					args.removeIf(x->x.startsWith("rankingfeatures=") || x.startsWith("scalarscale=") ||
						x.startsWith("scalarbounds=") || x.startsWith("scalarboundssha80=") || x.startsWith("cosine="));
					final boolean gated=mode.equals("full") || mode.equals("off");
					args.add("cosine="+(gated ? mode : "ranking"));
					args.add("rankingfeatures="+(gated ? "dimer" : mode));
				}
				boolean rejected=false;
				try{
					FamilyTailCascade.run(args.toArray(new String[0]));
				}catch(IllegalArgumentException failure){
					rejected=true;
				}
				check(rejected==(compact==1), "Unsupported compact mode must remain valid for legacy storage: "+mode);
				if(compact==1){
					check(!Files.exists(out) && !Files.exists(audit) && !Files.exists(rank),
						"Compact mode rejection precedes every writer: "+mode);
				}
			}
		}
	}
	/** Build matched three-query cascade arguments for 400-dimer versus 20-monomer storage. */
	private static String[] compactSearchArgs(Path reference, Path query, Path out, Path audit, Path bounds, int threads){
		return new String[]{"familyreferences="+reference, "familyreferencessha80="+hash(reference), "in="+query,
			"insha80="+hash(query),
			"out="+out, "auditout="+audit, "maxreads=3", "top=2", "k=5", "seed=0", "width=1024", "ddlreject=1", "ddlaccept=1025",
			"positional=f", "scalargates=t", "scalarz=3.8906", "cosine=ranking", "rankingfeatures=combined25", "scalarscale=0.1",
			"scalarbounds="+bounds, "scalarboundssha80="+hash(bounds), "t="+threads};
	}

	/*--------------------------------------------------------------*/
	/*---------------- Upper evidence and overrides ----------------*/
	/*--------------------------------------------------------------*/

	/** Additional negative samples widen raw evidence; a new lower may loosen a formerly lower-bound final gate. */
	private static void testAdditional(Path root) throws Exception{
		final Fixture f=new Fixture(root.resolve("additional"));
		cell(f.lower.resolve("cutoffs.tsv"), 1, 8, "0.6");
		f.bind();
		final Path old=f.root.resolve("old");
		FamilyReferenceExport.run(f.args(old));
		check(firstFamily(old).blosumFinal>=0.6, "Initial fixture final is bound by the old member lower");
		cell(f.lower.resolve("cutoffs.tsv"), 1, 8, "0.05");
		f.bind();
		cell(f.refIndex, 1, 10, "NEW_RETAINED_SCORES");
		final Path supplement=f.root.resolve("old_union.tsv"), extra=f.root.resolve("new_sample.tsv");
		write(supplement, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n0\tf0\t3\t90\t0.1\n");
		for(int observed=0; observed<2; observed++){
			write(extra, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n"+(observed==0 ? "" : "0\tf0\t500\t85\t0.2\n"));
			final Path out=f.root.resolve("new_"+observed);
			final java.util.ArrayList<String> args=new java.util.ArrayList<String>(
				Arrays.asList(supplementArgs(f, out, supplement)));
			args.add("additionalupper="+extra);
			args.add("additionaluppersha80="+hash(extra));
			args.add("additionalqueries=1000000");
			args.add("blosumfloor=0.130");
			FamilyReferenceExport.run(args.toArray(new String[0]));
			final FamilyReference actual=firstFamily(out);
			final double expected=observed==0 ? 0.130 : 1.1*0.2;
			check(actual.identityUpper>=1.1*90 && (double)Math.nextDown(actual.identityUpper)<1.1*90,
				"Missing or weaker new identity evidence preserves the approved old union maximum");
			check(actual.blosumFinal>=expected && (double)Math.nextDown(actual.blosumFinal)<expected &&
				actual.blosumFinal<firstFamily(old).blosumFinal,
				"Fresh lower replaces old binding lower while raw maxima union and global floor are applied once");
			check(actual.blosumLower<=0.05 && (double)Math.nextUp(actual.blosumLower)>0.05,
				"New member lower is stored with directed lower rounding");
			check(read(out.resolve("recipe.txt")).contains("additional_queries=1000000\n") &&
				read(out.resolve("recipe.txt")).contains("additional_families="+observed+"\n"),
				"No-observation families are distinguished from additional observed maxima");
			final List<String> quantiles=Files.readAllLines(out.resolve("gate_quantiles.tsv"), StandardCharsets.US_ASCII);
			check(quantiles.size()==7, "Exporter summarizes six raw/final metric distributions");
			final double[] expectedValues={30, 0.05, 90, observed==0 ? 0.1 : 0.2, actual.identityUpper, actual.blosumFinal};
			for(int row=0; row<6; row++){
				final String[] fields=quantiles.get(row+1).split("\t");
				check(fields[1].equals("1"), "Empty roster rows excluded from gate distributions");
				for(int col=2; col<fields.length; col++){
					check(Double.parseDouble(fields[col])==expectedValues[row],
						"Singleton quantiles preserve raw/final distinction");
				}
			}
		}
		final String[] negativeFloor=Arrays.copyOf(f.args(f.root.resolve("bad_floor")), 9);
		negativeFloor[8]="blosumfloor=-0.1";
		reject(negativeFloor);
		testAdditionalList(f, supplement, extra);
	}
	/** Verify per-metric maxima across listed upper tables and reject malformed/duplicate lists. */
	private static void testAdditionalList(Fixture f, Path first, Path second) throws Exception{
		write(second, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n0\tf0\t4\t95\t0.05\n");
		final Path list=f.root.resolve("additional_list.tsv");
		final String rowA=first+"\t"+hash(first)+"\t1000000\n", rowB=second+"\t"+hash(second)+"\t1000000\n";
		write(list, "path\tsha80\tnegative_count\n"+rowA+rowB);
		final Path out=f.root.resolve("listed");
		FamilyReferenceExport.run(listArgs(f, out, list));
		final FamilyReference actual=firstFamily(out);
		check(actual.identityUpper>=1.1*95 && (double)Math.nextDown(actual.identityUpper)<1.1*95,
			"Listed sources retain the stronger identity evidence from the second table");
		check(actual.blosumFinal>=1.1*.1 && (double)Math.nextDown(actual.blosumFinal)<1.1*.1,
			"Listed sources retain stronger BLOSUM evidence from the first table");
		final String recipe=read(out.resolve("recipe.txt"));
		check(recipe.contains("additional_tables=2\n") && recipe.contains("additional_families=1\n") &&
			recipe.contains("additional_pair_observations=7\n"),
			"Distinct families and per-source observations are counted separately; no unique pair claim");
		for(String row:new String[]{rowA+rowA, first+"\t"+hash(first)+"\t0\n", first+"\t00000000000000000000\t1000000\n", ""}){
			write(list, "path\tsha80\tnegative_count\n"+row);
			reject(listArgs(f, f.root.resolve("bad_list"), list));
		}
	}
	/** Extend fresh fixture arguments with the pinned additional-upper manifest. */
	private static String[] listArgs(Fixture f, Path out, Path list){
		final String[] base=f.args(out), args=Arrays.copyOf(base, base.length+2);
		args[base.length]="additionalupperlist="+list;
		args[base.length+1]="additionalupperlistsha80="+hash(list);
		return args;
	}

	/** Two round outputs must retain all prior increases and every non-acceptance payload byte. */
	private static void testRaise(Fixture fixture, Path baseline) throws Exception{
		final FamilyReference f=firstFamily(baseline);
		final double[] bounds={-2, 2, -3, 3, -4, 4, -5, 5, -6, 6};
		final FamilyReference scalar=new FamilyReference(f.id, f.rank, 30, f.residues, f.validWords, f.consensus, f.ddl,
			f.spectrum,
			f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept, f.identityUpper, f.blosumFinal,
			f.consensusRaw, 31, bounds);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA, 5, 1024, 0, 2, 1, 1,
			FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906);
		final Path parent=fixture.root.resolve("raise_parent.tsv"), table=fixture.root.resolve("raise_overrides.tsv");
		final ByteBuilder bytes=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		FamilyReferenceCodec.appendFiveLines(scalar, bytes);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty", 1), bytes);
		write(parent, bytes.toString());
		write(table, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n");
		final Path noChange=fixture.root.resolve("raise_empty");
		FamilyReferenceExport.run(raiseArgs(parent, table, noChange));
		check(Arrays.equals(Files.readAllBytes(parent), Files.readAllBytes(noChange.resolve("references.tsv"))),
			"Empty override table preserves every reference byte");
		write(table, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n0\tf0\t88\t0.130000001\n");
		final Path round1=fixture.root.resolve("raise_one");
		FamilyReferenceExport.run(raiseArgs(parent, table, round1));
		final FamilyReference after=firstFamily(round1);
		check(after.identityUpper==88 && after.blosumFinal>=0.130000001 && (double)Math.nextDown(after.blosumFinal)<0.130000001,
			"Raising BLOSUM uses directed rounding once and preserves identity");
		final List<String> oldLines=Files.readAllLines(parent, StandardCharsets.US_ASCII),
			newLines=Files.readAllLines(round1.resolve("references.tsv"), StandardCharsets.US_ASCII);
		check(oldLines.size()==newLines.size(), "Raise preserves physical five-line record census");
		for(int i=0; i<oldLines.size(); i++){
			if(i!=6){
				check(oldLines.get(i).equals(newLines.get(i)), "Raise preserves header/metadata/payload line "+i);
			}else{
				final String[] a=oldLines.get(i).split("\t"), b=newLines.get(i).split("\t");
				for(int j=0; j<a.length; j++){
					if(j!=5){check(a[j].equals(b[j]), "Raise preserves non-BLOSUM gate/scalar field "+j);}
				}
			}
		}
		check(read(round1.resolve("summary.tsv")).startsWith("rows_read\t") &&
			read(round1.resolve("summary.tsv")).contains("\n1\t0\t1\t1\t0\t"),
			"One BLOSUM-only change is reported separately from identity changes");
		write(table, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n0\tf0\t89\t"+Double.toString((double)after.blosumFinal)+"\n");
		final Path round2=fixture.root.resolve("raise_two");
		FamilyReferenceExport.run(raiseArgs(round1.resolve("references.tsv"), table, round2));
		check(firstFamily(round2).identityUpper==89 && firstFamily(round2).blosumFinal==after.blosumFinal,
			"Next round preserves prior BLOSUM increase while raising identity");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(round2.resolve("references.tsv").toString(), 1, 2);
		check(loaded.header.scalarZ==3.8906 && loaded.all[0].scalarMask==31 && Arrays.equals(loaded.all[0].scalarBounds, bounds),
			"Published raised package reloads with unchanged enabled scalar gates");
		int index=0;
		for(String row:new String[]{"0\tf0\t87\t0.2", "0\tf0\t88\t0", "1\tempty\t88\t0.2", "0\twrong\t88\t0.2",
			"2\tf0\t88\t0.2", "0\tf0\t1e309\t0.2", "0\tf0\t88\t1e309", "0\tf0\t88\tNA",
			"0\tf0\t111\t0.2", "0\tf0\t88\t0.2\n0\tf0\t88\t0.3"}){
			write(table, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n"+row+"\n");
			final Path bad=fixture.root.resolve("raise_bad_"+(index++));
			reject(raiseArgs(parent, table, bad));
			check(!Files.exists(bad), "Invalid override must not publish output");
		}
		write(table, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n");
		final String[] wrongPin=raiseArgs(parent, table, fixture.root.resolve("raise_bad_pin"));
		replace(wrongPin, "parentsha80", "00000000000000000000");
		reject(wrongPin);
		reject(raiseArgs(parent, table, round1));
		check(Files.exists(round1.resolve("PASS")), "Existing round preserved");
		testCandidateExports(fixture.root, parent, scalar);
	}

	/** Historical references are final floors only; current export owns every other payload field. */
	private static void testHistoricalFloors(Fixture fixture, Path baseline) throws Exception{
		final FamilyReference before=firstFamily(baseline);
		final Path full=fixture.root.resolve("historical_full.tsv"), compact=fixture.root.resolve("historical_compact.tsv");
		writeFloor(full, false, new FamilyReference(before.id, before.rank, 9, before.residues, before.validWords,
			before.consensus,
			before.ddl, before.spectrum, before.cosineLower, 1, 0, before.cosineAccept, 89, 0, before.consensusRaw),
			FamilyReference.empty("empty", 1));
		writeFloor(compact, true, new FamilyReference(before.id, before.rank, 9, before.residues, before.validWords,
			before.consensus,
			before.ddl, before.spectrum, before.cosineLower, 1, 0, before.cosineAccept, before.identityUpper, 0.2f, before.consensusRaw),
			FamilyReference.empty("empty", 1));
		for(int variant=0; variant<2; variant++){
			final Path floor=variant==0 ? full : compact, out=fixture.root.resolve("historical_floor_"+variant);
			FamilyReferenceExport.run(floorArgs(fixture, out, floor));
			final FamilyReference actual=firstFamily(out);
			check(actual.identityUpper==(variant==0 ? 89 : before.identityUpper) &&
				actual.blosumFinal==(variant==0 ? before.blosumFinal : 0.2f),
					"Historical final floors max independently by field");
			check(actual.identityLower==before.identityLower && actual.blosumLower==before.blosumLower &&
				actual.cosineAccept==before.cosineAccept &&
				Arrays.equals(actual.spectrum, before.spectrum) && Arrays.equals(actual.ddl, before.ddl) &&
					Arrays.equals(actual.consensusRaw, before.consensusRaw),
				"Historical floors preserve every current non-final-gate payload");
			check(read(out.resolve("recipe.txt")).contains("historical_floor_reference="+floor+"\n") &&
				read(out.resolve("inputs.tsv")).contains(floor+"\t"+hash(floor)),
					"Historical floor input and policy are recorded");
		}
		final Path equal=fixture.root.resolve("historical_equal.tsv"), equalOut=fixture.root.resolve("historical_equal");
		writeFloor(equal, false, before, FamilyReference.empty("empty", 1));
		FamilyReferenceExport.run(floorArgs(fixture, equalOut, equal));
		check(Arrays.equals(Files.readAllBytes(baseline.resolve("references.tsv")), Files.readAllBytes(equalOut.resolve("references.tsv"))),
			"Equal historical floors preserve emitted reference bytes");
		int index=0;
		final byte[] raw=before.consensusRaw;
		for(FamilyReference[] rows:new FamilyReference[][]{
			{new FamilyReference("wrong", 0, 9, before.residues, before.validWords, new ProteinSequence("wrong", raw),
				before.ddl, before.spectrum,
				before.cosineLower, 1, 0, before.cosineAccept, 89, 0, raw), FamilyReference.empty("empty", 1)},
			{
				before, new FamilyReference("empty", 1, 1, 1, 1, new ProteinSequence("empty", raw), before.ddl, before.spectrum,
				before.cosineLower, 1, 0, before.cosineAccept, 89, 0, raw)
			},
			{
				new FamilyReference(before.id, 0, 9, before.residues, before.validWords,
					new ProteinSequence(before.id, "CCCC".getBytes(StandardCharsets.US_ASCII)),
				before.ddl, before.spectrum, before.cosineLower, 1, 0, before.cosineAccept, 89, 0, "CCCC".getBytes(StandardCharsets.US_ASCII)),
				FamilyReference.empty("empty", 1)
			}}){
			final Path bad=fixture.root.resolve("historical_bad_"+index+".tsv"),
				badOut=fixture.root.resolve("historical_bad_"+(index++));
			writeFloor(bad, false, rows);
			reject(floorArgs(fixture, badOut, bad));
			check(!Files.exists(badOut), "Mismatched historical floor must not publish output");
		}
		final String[] badPin=floorArgs(fixture, fixture.root.resolve("historical_bad_pin"), equal);
		replace(badPin, "floorreferencesha80", "00000000000000000000");
		reject(badPin);
	}

	/** Extend fresh fixture arguments with one pinned historical final-floor reference. */
	private static String[] floorArgs(Fixture f, Path out, Path floor){
		final String[] base=f.args(out), args=Arrays.copyOf(base, base.length+2);
		args[base.length]="floorreference="+floor;
		args[base.length+1]="floorreferencesha80="+hash(floor);
		return args;
	}

	/** Serialize borrowed records into a synthetic legacy or projected-monomer floor file.
	 * Projection creates a new live record and spectrum; borrowed input records remain unchanged. */
	private static void writeFloor(Path path, boolean compact, FamilyReference... records) throws Exception{
		int live=0, empty=0;
		for(FamilyReference record:records){
			if(record.members>0){live++;}else{empty++;}
		}
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(compact ? FamilyReferenceCodec.MONOMER_SCHEMA : FamilyReferenceCodec.SCHEMA,
			5, 1024, 0, records.length, live, empty, FamilyReferenceCodec.FORMULA_VERSION,
				FamilyReferenceCodec.SPECTRUM_ENCODING, 0);
		final ByteBuilder text=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		for(FamilyReference record:records){
			FamilyReference out=record;
			if(compact && record.members>0){
				final float[] projected=new float[FamilyReferenceCodec.MONOMER_WIDTH];
				FamilyTailCascade.dimerRowSums(record.spectrum, projected);
				out=new FamilyReference(record.id, record.rank, record.members, record.residues, record.validWords,
					record.consensus,
					record.ddl, projected, record.cosineLower, record.identityLower, record.blosumLower, record.cosineAccept,
					record.identityUpper, record.blosumFinal, record.consensusRaw, record.scalarMask, record.scalarBounds);
			}
			FamilyReferenceCodec.appendFiveLines(out, text);
		}
		write(path, text.toString());
	}
	/** Multi-seed tables share one parent; global-floor candidates change only acceptance bounds. */
	private static void testCandidateExports(Path root, Path parent, FamilyReference before) throws Exception{
		final Path a=root.resolve("candidate_a.tsv"), b=root.resolve("candidate_b.tsv"),
			empty=root.resolve("candidate_empty.tsv"), list=root.resolve("candidate_list.tsv");
		write(a, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n0\tf0\t89\t"+(double)before.blosumFinal+"\n");
		write(b, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n0\tf0\t88\t0.200000001\n");
		write(empty, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n");
		final String rowA=a+"\t"+hash(a)+"\n", rowB=b+"\t"+hash(b)+"\n", rowEmpty=empty+"\t"+hash(empty)+"\n";
		Path prior=null;
		for(int order=0; order<2; order++){
			write(list, "path\tsha80\n"+(order==0 ? rowA+rowB+rowEmpty : rowEmpty+rowB+rowA));
			final Path out=root.resolve("candidate_union_"+order);
			FamilyReferenceExport.run(candidateArgs(parent, out, "overrideslist="+list, "overrideslistsha80="+hash(list)));
			final FamilyReference f=firstFamily(out);
			check(f.identityUpper==89 && f.blosumFinal>=0.200000001 && (double)Math.nextDown(f.blosumFinal)<0.200000001,
				"Opposing per-seed changes merge against one parent, not sequential intermediate gates");
			check(read(out.resolve("recipe.txt")).contains("override_tables=3\nsource_rows=2\n") &&
				read(out.resolve("summary.tsv")).contains("\n1\t1\t1\t1\t0\t"),
					"Counts separate source rows from unique changed families");
			check(read(out.resolve("merged_overrides.tsv")).contains("0\tf0\t89.0\t0.200000001\n"),
				"Merged raw-double override is retained before final rounding");
			check(read(out.resolve("override_sources.tsv")).contains("0\tf0\t"+a+"\t"+b+"\n"),
				"Each merged field names the source that supplied its maximum");
			if(prior!=null){
				check(Arrays.equals(Files.readAllBytes(prior), Files.readAllBytes(out.resolve("references.tsv"))),
					"Override union is independent of table order");
			}
			prior=out.resolve("references.tsv");
		}
		final Path floors=root.resolve("candidate_floors");
		FamilyReferenceExport.run(candidateArgs(parent, floors, "identityfloor=40", "blosumfloor=0.16"));
		final FamilyReference f=firstFamily(floors);
		check(f.identityUpper==before.identityUpper && f.blosumFinal>=.16 && (double)Math.nextDown(f.blosumFinal)<.16,
			"Global floors preserve stricter parent gates and upward-round a raising floor");
		check(f.identityLower==before.identityLower && f.blosumLower==before.blosumLower && f.scalarMask==before.scalarMask &&
			Arrays.equals(f.scalarBounds, before.scalarBounds) && Arrays.equals(f.consensusRaw, before.consensusRaw) &&
			Arrays.equals(f.ddl, before.ddl) && Arrays.equals(f.spectrum, before.spectrum),
				"Floor candidate preserves every biological/scalar payload");
		final Path both=root.resolve("candidate_both_floors");
		FamilyReferenceExport.run(candidateArgs(parent, both, "identityfloor=90", "blosumfloor=0.16"));
		check(firstFamily(both).identityUpper==90, "Identity max-floor can independently raise a parent gate");
		write(b, FamilyReferenceExport.ACCEPT_OVERRIDE_HEADER+"\n0\tf0\t87\t0.3\n");
		write(list, "path\tsha80\n"+rowA+b+"\t"+hash(b)+"\n");
		reject(candidateArgs(parent, root.resolve("candidate_bad_parent"), "overrideslist="+list,
			"overrideslistsha80="+hash(list)));
		for(String text:new String[]{"path\tsha80\n", "path\tsha80\n"+rowA+rowA, "path\tsha80\n"+a+"\t00000000000000000000\n"}){
			write(list, text);
			reject(candidateArgs(parent, root.resolve("candidate_bad_list"), "overrideslist="+list,
				"overrideslistsha80="+hash(list)));
		}
		reject(candidateArgs(parent, root.resolve("candidate_bad_floor"), "identityfloor=101"));
		reject(candidateArgs(parent, root.resolve("candidate_bad_blosum"), "blosumfloor=-0.1"));
	}
	/** Build fresh raise-mode arguments and append caller-provided floor or manifest options. */
	private static String[] candidateArgs(Path parent, Path out, String... extra){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>(Arrays.asList("mode=raise", "parent="+parent,
			"parentsha80="+hash(parent), "out="+out, "maxfamilies=10", "threads=1"));
		args.addAll(Arrays.asList(extra));
		return args.toArray(new String[0]);
	}
	/** Build fresh raise-mode arguments bound to one parent and one override table. */
	private static String[] raiseArgs(Path parent, Path table, Path out){
		return new String[]{"mode=raise", "parent="+parent, "parentsha80="+hash(parent), "overrides="+table,
			"overridessha80="+hash(table), "out="+out, "maxfamilies=10", "threads=1"};
	}

	/*--------------------------------------------------------------*/
	/*---------------- Scalar and prepared records  ----------------*/
	/*--------------------------------------------------------------*/

	/** Check scalar-table bounds, disabled zero-SD/small-family gates and malformed rows.
	 * The final one-member CLI case checks header binding without enabling scalar gates. */
	private static void testScalars(Fixture fixture, Path baseline) throws Exception{
		final FamilyReference f=firstFamily(baseline);
		final StringBuilder header=new StringBuilder("family_rank\tfamily_id\tN");
		for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
			for(String metric:new String[]{"mean", "sd", "min", "max"}){
				header.append('\t').append(name).append('_').append(metric);
			}
		}
		final StringBuilder text=new StringBuilder(header).append("\n0\tf0\t30");
		for(int i=0; i<5; i++){text.append(i==1 ? "\t2\t0\t2\t2" : "\t2\t1\t0\t4");}
		text.append("\n1\tempty\t0");
		for(int i=0; i<20; i++){text.append("\tNA");}
		text.append('\n');
		final Path table=fixture.root.resolve("scalars.tsv");
		write(table, text.toString());
		for(double z:new double[]{3.2905, 3.8906}){
			final FamilyReference[] records={new FamilyReference(f.id, 0, 30, f.residues, f.validWords, f.consensus,
				f.ddl, f.spectrum,
				f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept, f.identityUpper, f.blosumFinal, f.consensusRaw), FamilyReference.empty("empty", 1)};
			FamilyReferenceExport.applyScalarBounds(table, records, z);
			check(records[0].scalarMask==29 && records[1].scalarMask==0, "Zero SD and N0 disable gates without nonfinite values");
			check(Math.abs(records[0].scalarBounds[0]-(z==3.2905 ? -1.2905 : -1.8906))<1e-14 &&
				Math.abs(records[0].scalarBounds[1]-(z==3.2905 ? 5.2905 : 5.8906))<1e-14,
					"Both population-SD bounds match independent expected values");
			check(records[0].identityUpper==f.identityUpper && records[0].blosumFinal==f.blosumFinal && records[0].ddl==f.ddl,
				"Scalar export preserves alignment gates and ranking arrays");
		}
		for(String mutation:new String[]{"small", "mismatch", "duplicate", "missing", "nonfinite", "negative_sd"}){
			final long n=mutation.equals("small") ? 29 : 30;
			final FamilyReference[] records={new FamilyReference(f.id, 0, n, f.residues, f.validWords, f.consensus,
				f.ddl, f.spectrum,
				f.cosineLower, f.identityLower, f.blosumLower, f.cosineAccept, f.identityUpper, f.blosumFinal, f.consensusRaw), FamilyReference.empty("empty", 1)};
			String changed=text.toString();
			if(mutation.equals("small")){changed=changed.replace("\tf0\t30", "\tf0\t29");}
			if(mutation.equals("mismatch")){changed=changed.replace("\tf0\t30", "\tf0\t31");}
			if(mutation.equals("duplicate")){changed=changed.replace("1\tempty", "0\tempty");}
			if(mutation.equals("missing")){changed=changed.substring(0, changed.indexOf("1\tempty"));}
			if(mutation.equals("nonfinite")){changed=changed.replace("\t2\t1\t0\t4", "\t1e309\t1\t0\t4");}
			if(mutation.equals("negative_sd")){changed=changed.replace("\t2\t1\t0\t4", "\t2\t-1\t0\t4");}
			write(table, changed);
			boolean failed=false;
			try{
				FamilyReferenceExport.applyScalarBounds(table, records, 3.2905);
			}catch(IllegalArgumentException expected){
				failed=true;
			}
			check(mutation.equals("small") ? !failed && records[0].scalarMask==0 : failed, "Scalar table case: "+mutation);
		}
		// Full CLI wiring uses the existing one-member fixture: all five gates must be disabled.
		write(table, text.toString().replace("\tf0\t30", "\tf0\t1").replace("\t2\t1\t0\t4", "\t2\t0\t2\t2"));
		final String[] base=fixture.args(fixture.root.resolve("scalar_export")), args=Arrays.copyOf(base, base.length+3);
		args[base.length]="scalars="+table;
		args[base.length+1]="scalarssha80="+hash(table);
		args[base.length+2]="scalarz=3.2905";
		FamilyReferenceExport.run(args);
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(fixture.root.resolve("scalar_export/references.tsv").toString(), 1, 2);
		check(loaded.header.scalarZ==3.2905 && loaded.all[0].scalarMask==0, "CLI scalar pin/z reaches packed reference header");
	}

	/** Rounded moments reproduce the table path bit-for-bit, including the zero-SD and minimum-count boundaries. */
	private static void testScalarApi(Fixture fixture, Path baseline) throws Exception{
		final FamilyReference base=firstFamily(baseline);
		final double[] expected={0.123456789012, 0.111111111111, -0.222222222222, 0.333333333333};
		for(int members:new int[]{1, 29, 30}){
			final FamilyMemberScalars.Moments moments=new FamilyMemberScalars.Moments(1);
			moments.n[0]=members;
			final ByteBuilder table=new ByteBuilder("family_rank\tfamily_id\tN");
			for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
				for(String metric:new String[]{"mean", "sd", "min", "max"}){
					table.tab().append(name).append('_').append(metric);
				}
			}
			table.nl().append("0\tf0\t").append(members);
			for(int i=0; i<5; i++){
				moments.mean[i]=0.12345678901234;
				moments.min[i]=-0.22222222222224;
				moments.max[i]=0.33333333333334;
				final double sd=members==1 ? 0 : i==1 ? 1e-14 : 0.11111111111114;
				moments.m2[i]=members*sd*sd;
				table.tab().append(Double.toString(expected[0])).tab()
					.append(Double.toString(members==1 || i==1 ? 0 : expected[1]))
					.tab().append(Double.toString(expected[2])).tab().append(Double.toString(expected[3]));
			}
			final Path path=fixture.root.resolve("scalar_api_"+members+".tsv");
			write(path, table.nl().toString());
			final FamilyReference original=new FamilyReference(base.id, 0, members, base.residues, base.validWords,
				base.consensus,
				base.ddl, base.spectrum, base.cosineLower, base.identityLower, base.blosumLower, base.cosineAccept,
					base.identityUpper, base.blosumFinal, base.consensusRaw);
			final double[] statistics=new double[20];
			moments.exportStatistics(0, statistics, new ByteBuilder());
			for(double z:new double[]{3.2905, 3.8906}){
				final FamilyReference[] fromFile={original};
				FamilyReferenceExport.applyScalarBounds(path, fromFile, z);
				final FamilyReference direct=FamilyReferenceExport.withScalarBounds(original, statistics, z);
				check(Arrays.equals(FamilyReferenceCodec.appendFiveLines(fromFile[0], new ByteBuilder()).toBytes(),
					FamilyReferenceCodec.appendFiveLines(direct, new ByteBuilder()).toBytes()),
						"File and in-memory scalar gates preserve every encoded byte");
				check(direct.scalarMask==(members<30 ? 0 : 29), "Small families and decimal-rounded zero SD remain disabled");
				if(members==30){
					check(Double.doubleToLongBits(direct.scalarBounds[0])==Double.doubleToLongBits(expected[0]-z*expected[1]) &&
						Double.doubleToLongBits(direct.scalarBounds[1])==Double.doubleToLongBits(expected[0]+z*expected[1]), "Bounds use the literal rounded statistics");
					check(direct.scalarBounds[0]!=moments.mean[0]-z*Math.sqrt(moments.m2[0]/members),
						"Fixture distinguishes unrounded in-memory bounds");
				}
				check(original.scalarMask==0 && direct.ddl==original.ddl && direct.consensus==original.consensus,
					"Binding scalar gates leaves the input and other payloads untouched");
			}
		}
		final FamilyReference empty=FamilyReference.empty("empty", 0);
		check(FamilyReferenceExport.withScalarBounds(empty, null, 3.8906)==empty,
			"Empty family remains the same explicit empty record");
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(empty, new double[20], 3.8906));
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, null, 3.8906));
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, new double[19], 3.8906));
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, new double[20], 0));
		for(double bad:new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}){
			final double[] invalid=new double[20];
			invalid[0]=bad;
			rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, invalid, 3.8906));
		}
		final double[] invalid=new double[20];
		invalid[1]=1;
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, invalid, 3.8906));
		invalid[1]=-1;
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(base, invalid, 3.8906));
		final FamilyReference live=new FamilyReference(base.id, 0, 30, base.residues, base.validWords, base.consensus,
			base.ddl, base.spectrum, base.cosineLower, base.identityLower, base.blosumLower, base.cosineAccept,
				base.identityUpper, base.blosumFinal, base.consensusRaw);
		final double[] extremes=new double[20];
		extremes[0]=extremes[2]=extremes[3]=1e300;
		extremes[1]=1;
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(live, extremes, 3.8906));
		extremes[0]=extremes[2]=extremes[3]=Double.MAX_VALUE;
		extremes[1]=Double.MAX_VALUE;
		rejectScalar(()->FamilyReferenceExport.withScalarBounds(live, extremes, 3.8906));
	}
	/** The existing file exporter is the oracle for complete records, including scalar/rounding boundaries. */
	private static void testPreparedReferences(Path root) throws Exception{
		final FamilyMemberPreparation worker=new FamilyMemberPreparation(4, 2, null,
			FamilyCalibration.MemberRetention.METRIC_SPECIFIC);
		final byte[] consensus=new byte[100];
		Arrays.fill(consensus, (byte)'A');
		final byte[] alphabet="ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII);
		final int[] sizes={1, 30, 400};
		final double[][] maxima={{.5, 80, -.001, 0, 0},{.99, 70, .12, .13, 3.8906},{1, 100, 2, 0, 3.8906}};
		FamilyMemberPreparation.Prepared live=null;
		for(int test=0; test<sizes.length; test++){
			final int n=sizes[test];
			final byte[][] members=new byte[n][];
			final DDLFamilyProfiles.Profile votes=new DDLFamilyProfiles.Profile(1024);
			final char[] sketch=new char[1024];
			for(int i=0; i<n; i++){
				members[i]=consensus.clone();
				members[i][99]=alphabet[i%alphabet.length];
				DDLFamilyProfiles.sketch(members[i], 5, 0, sketch);
				votes.add(sketch);
			}
			final FamilyMemberPreparation.Prepared p=worker.prepare(0, "f0", members, consensus);
			live=p;
			final Fixture f=new Fixture(root.resolve("prepared_reference_"+n));
			write(f.refs.resolve("members.tsv"),
				"roster_rank\tfamily_id\tmembers\tresidues\tvalid_words\n0\tf0\t"+n+"\t"+p.residues+"\t"+p.validWords+"\n1\tempty\t0\t0\t0\n");
			final ByteBuilder profiles=new ByteBuilder(DDLFamilyProfiles.HEADER16+"\n");
			for(boolean literal:new boolean[]{true, false}){
				profiles.append("tail\tf0\t1024\t").append(n).tab().append(p.residues).tab().append(p.validWords).tab()
					.append(literal ? "literal" : "nonempty");
				for(int field=0; field<3; field++){
					profiles.tab();
					for(int b=0; b<1024; b++){
						if(b>0){profiles.append(',');}
						profiles.append(field==0 ? (int)votes.value(b, literal) :
							field==1 ? votes.votes(b, literal) : votes.occupied(b));
					}
				}
				profiles.nl();
			}
			write(f.refs.resolve("profiles.tsv"), profiles.toString());
			final ByteBuilder spectrum=new ByteBuilder();
			ProteinSpectrum.writeRow(spectrum, "tail", "f0", p.spectrum);
			write(f.refs.resolve("spectra.tsv"), ProteinSpectrum.header(2)+"\n"+spectrum);
			cell(f.refs.resolve("summary.tsv"), 1, 9, Integer.toString(n));
			cell(f.refs.resolve("summary.tsv"), 1, 10, Long.toString(p.residues));
			final ByteBuilder lower=new ByteBuilder(FamilyCalibration.CUTOFF_HEADER+"\n0\tf0\t").append(n).append("\tCALIBRATED");
			for(FamilyCalibration.Quantile q:p.lower){
				lower.tab().append(Double.toString(q.cutoff)).tab().append(q.retained);
			}
			write(f.lower.resolve("cutoffs.tsv"), lower.append("\n1\tempty\t0\tNO_MEMBERS\tNA\t0\tNA\t0\tNA\t0\n").toString());
			final double[] upper=maxima[test];
			for(int metric=0; metric<3; metric++){
				cell(f.upper.resolve("uppers.tsv"), 1, metric==0 ? 3 : metric+4, Double.toString(upper[metric]));
			}
			f.bind(n, p.residues);
			final Path target=f.root.resolve("file_export"), scalars=f.root.resolve("scalars.tsv");
			final ByteBuilder table=new ByteBuilder("family_rank\tfamily_id\tN");
			for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
				for(String metric:new String[]{"mean", "sd", "min", "max"}){
					table.tab().append(name).append('_').append(metric);
				}
			}
			table.nl().append("0\tf0\t").append(n);
			for(double statistic:p.statistics){table.tab().append(Double.toString(statistic));}
			table.append("\n1\tempty\t0");
			for(int i=0; i<20; i++){table.append("\tNA");}
			table.nl();
			write(scalars, table.toString());
			final String[] base=f.args(target), args=Arrays.copyOf(base, base.length+(upper[4]>0 ? 4 : 1));
			args[base.length]="blosumfloor="+upper[3];
			if(upper[4]>0){
				args[base.length+1]="scalars="+scalars;
				args[base.length+2]="scalarssha80="+hash(scalars);
				args[base.length+3]="scalarz="+upper[4];
			}
			FamilyReferenceExport.run(args);
			final FamilyReference direct=FamilyReferenceExport.fromPrepared(p, upper[0], upper[1], upper[2], upper[3], upper[4]);
			final List<String> rows=Files.readAllLines(target.resolve("references.tsv"), StandardCharsets.US_ASCII);
			final byte[] expected=(String.join("\n", rows.subList(5, 10))+"\n").getBytes(StandardCharsets.US_ASCII);
			check(Arrays.equals(expected, FamilyReferenceCodec.appendFiveLines(direct, new ByteBuilder()).toBytes()),
				"Prepared final reference must match every file-export byte at member count"+n);
			check(direct.scalarMask==firstFamily(target).scalarMask && (n==1 || direct.scalarMask!=0),
				"Scalar enablement must match the requested measured population");
		}
		final FamilyMemberPreparation.Prepared family=live, empty=worker.prepare(1, "empty", new byte[0][], consensus);
		rejectScalar(()->FamilyReferenceExport.fromPrepared(null, 0, 0, 0, 0, 0));
		rejectScalar(()->FamilyReferenceExport.fromPrepared(empty, 0, 0, 0, 0, 0));
		for(double bad:new double[]{Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY}){
			rejectScalar(()->FamilyReferenceExport.fromPrepared(family, bad, 0, 0, 0, 0));
			rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, bad, 0, 0, 0));
			rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 0, bad, 0, 0));
			rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 0, 0, bad, 0));
			rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 0, 0, 0, bad));
		}
		rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 1.01, 0, 0, 0, 0));
		rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 101, 0, 0, 0));
		rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 0, 0, -1, 0));
		rejectScalar(()->FamilyReferenceExport.fromPrepared(family, 0, 0, 0, 0, -1));
		System.out.println("FAMILY_PREPARED_REFERENCE_PASS file_bytes=true scalar_bounds=true raw_rounding=true invalid_rejected=true");
	}

	/** Count an expected IllegalArgumentException; any other failure propagates.
	 * A normal return is a failed rejection test. */
	private static void rejectScalar(Runnable action){
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			checks++;
			return;
		}
		throw new AssertionError("Invalid scalar API input was accepted");
	}

	/*--------------------------------------------------------------*/
	/*---------------- Supplement and input checks  ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent raw-max union oracle: one metric can rise while the other retains its old bound. */
	private static void testSupplemental(Fixture f, Path baseline) throws Exception{
		final FamilyReference before=firstFamily(baseline);
		int index=0;
		for(double[] maxima:new double[][]{{40, 0.2},{90, -0.1},{80, -0.001}}){
			final Path supplement=f.root.resolve("supplement_"+index+".tsv"), out=f.root.resolve("union_"+(index++));
			write(supplement, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n0\tf0\t3\t"+maxima[0]+"\t"+maxima[1]+"\n");
			FamilyReferenceExport.run(supplementArgs(f, out, supplement));
			final FamilyReference actual=firstFamily(out);
			final double rawIdentity=Math.max(80, maxima[0]), rawBlosum=Math.max(-0.001, maxima[1]);
			final double id=Math.max(30, Math.max(1.1*rawIdentity, rawIdentity+2));
			final double bl=Math.max(0, Math.max(1.1*rawBlosum, rawBlosum+0.005));
			check(actual.identityUpper>=id && (double)Math.nextDown(actual.identityUpper)<id,
				"Union identity margin applied once to raw maximum");
			check(actual.blosumFinal>=bl && (double)Math.nextDown(actual.blosumFinal)<bl,
				"Union BLOSUM margin applied once to raw maximum");
			check(actual.identityLower==before.identityLower && actual.blosumLower==before.blosumLower &&
				actual.cosineAccept==before.cosineAccept,
				"Supplement changes only raw identity/BLOSUM uppers, never lowers or cosine");
			check(read(out.resolve("recipe.txt")).contains("supplemental_pair_observations=3\n") &&
				read(out.resolve("inputs.tsv")).contains(supplement+"\t"+hash(supplement)),
					"Supplement count and exact input pin retained");
		}
		final Path none=f.root.resolve("supplement_empty.tsv"), unchanged=f.root.resolve("union_empty");
		write(none, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n");
		FamilyReferenceExport.run(supplementArgs(f, unchanged, none));
		check(Arrays.equals(Files.readAllBytes(baseline.resolve("references.tsv")), Files.readAllBytes(unchanged.resolve("references.tsv"))),
			"No additional pairs preserves every packed reference byte");
		for(String row:new String[]{"2\tf0\t3\t90\t0.2", "1\tempty\t3\t90\t0.2", "0\twrong\t3\t90\t0.2",
			"0\tf0\t0\t90\t0.2", "0\tf0\t100001\t90\t0.2", "0\tf0\t3\t101\t0.2",
			"0\tf0\t3\tNaN\t0.2", "0\tf0\t3\t90\tInfinity", "0\tf0\t3\t90\t0.2\n0\tf0\t3\t90\t0.2"}){
			final Path supplement=f.root.resolve("bad_supplement_"+index+".tsv"), out=f.root.resolve("bad_union_"+(index++));
			write(supplement, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n"+row+"\n");
			reject(supplementArgs(f, out, supplement));
			check(!Files.exists(out), "Malformed supplement must not publish output");
		}
		final String[] badPin=supplementArgs(f, f.root.resolve("bad_supplement_pin"), none);
		replace(badPin, "supplementaluppersha80", "00000000000000000000");
		reject(badPin);
	}
	/** Extend fresh fixture arguments with one pinned sparse supplemental-upper table. */
	private static String[] supplementArgs(Fixture f, Path out, Path supplement){
		final String[] base=f.args(out), args=Arrays.copyOf(base, base.length+2);
		args[base.length]="supplementalupper="+supplement;
		args[base.length+1]="supplementaluppersha80="+hash(supplement);
		return args;
	}

	/** Complete packaged uppers replace bootstrap identity/BLOSUM, not cosine, and remain distinct from sparse unions. */
	private static void testAlignmentReplacement(Path root) throws Exception{
		final Fixture f=new Fixture(root.resolve("alignment_replacement"));
		cell(f.upper.resolve("uppers.tsv"), 1, 3, "0.90");
		cell(f.upper.resolve("uppers.tsv"), 1, 5, "99");
		cell(f.upper.resolve("uppers.tsv"), 1, 6, "1.0");
		f.bind();
		final Path table=f.root.resolve("alignment.tsv"), out=f.root.resolve("out");
		write(table, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n0\tf0\t500\t70\t0.04\n");
		FamilyReferenceExport.run(alignmentArgs(f, out, table, 1000000));
		final FamilyReference actual=firstFamily(out);
		check(actual.identityUpper>=77 && (double)Math.nextDown(actual.identityUpper)<77,
			"Complete alignment upper replaces a stronger bootstrap identity maximum instead of max-unioning it");
		check(actual.blosumFinal>=0.045 && (double)Math.nextDown(actual.blosumFinal)<0.045,
			"Complete alignment upper replaces a stronger bootstrap BLOSUM maximum before applying the unchanged formula");
		check(actual.cosineAccept>=0.99 && (double)Math.nextDown(actual.cosineAccept)<0.99,
			"Replacement leaves bootstrap cosine serialized for ranking-only provenance rather than fabricating cosine evidence");
		final String recipe=read(out.resolve("recipe.txt"));
		check(recipe.contains("alignment_queries=1000000\n") && recipe.contains("alignment_families=1\n") &&
			recipe.contains("alignment_pair_observations=500\n") &&
				recipe.contains("legacy_cosine_upper=retained_from_upperindex;not_newly_calibrated;requires_cosine_ranking_downstream\n"),
			"Replacement recipe records the complete panel and the inert cosine provenance");
		int index=0;
		for(String row:new String[]{"", "0\tf0\t500\t70\t0.04\n0\tf0\t500\t70\t0.04", "0\twrong\t500\t70\t0.04",
			"1\tempty\t500\t70\t0.04", "0\tf0\t0\t70\t0.04", "0\tf0\t1000001\t70\t0.04", "0\tf0\t500\tNaN\t0.04"}){
			final Path bad=f.root.resolve("bad_alignment_"+(index++)+".tsv"), badOut=f.root.resolve("bad_alignment_"+index);
			write(bad, FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n"+(row.isEmpty() ? "" : row+"\n"));
			reject(alignmentArgs(f, badOut, bad, 1000000));
			check(!Files.exists(badOut), "Malformed complete upper must not publish output");
		}
		reject(alignmentArgs(f, f.root.resolve("bad_alignment_population"), table, 0));
		final String[] badPin=alignmentArgs(f, f.root.resolve("bad_alignment_pin"), table, 1000000);
		replace(badPin, "alignmentuppersha80", "00000000000000000000");
		reject(badPin);
	}
	/** Verify one-million-query metadata is accepted without changing packed payloads,
	 * and reject mismatched population, alignment and top-count summaries. */
	private static void testMillionQueryBootstrap(Path root) throws Exception{
		final Fixture f=new Fixture(root.resolve("million_query_bootstrap"));
		final Path baseline=f.root.resolve("baseline");
		FamilyReferenceExport.run(f.args(baseline));
		cell(f.upper.resolve("uppers.tsv"), 1, 2, "1000000");
		f.writeUpperSummary(1000000, 500, 500);
		f.bind();
		final Path out=f.root.resolve("out");
		FamilyReferenceExport.run(f.args(out));
		check(read(out.resolve("PASS")).equals(FamilyReferenceExport.PASS+"\n"),
			"One-million upper populations must remain valid bootstrap inputs");
		check(Arrays.equals(Files.readAllBytes(baseline.resolve("references.tsv")), Files.readAllBytes(out.resolve("references.tsv"))),
			"Changing only the upper query population preserves every packed reference byte");
		final Fixture bad=new Fixture(root.resolve("million_query_mismatch"));
		cell(bad.upper.resolve("uppers.tsv"), 1, 2, "1000000");
		bad.bind();
		reject(bad.args(bad.root.resolve("out")));
		final Fixture badAligned=new Fixture(root.resolve("million_alignment_mismatch"));
		cell(badAligned.upper.resolve("uppers.tsv"), 1, 2, "1000000");
		badAligned.writeUpperSummary(1000000, 500, 499);
		badAligned.bind();
		reject(badAligned.args(badAligned.root.resolve("out")));
		final Fixture badTop=new Fixture(root.resolve("million_top_mismatch"));
		cell(badTop.upper.resolve("uppers.tsv"), 1, 2, "1000000");
		badTop.writeUpperSummary(1000000, 501, 500);
		badTop.bind();
		reject(badTop.args(badTop.root.resolve("out")));
	}
	/** Verify deterministic 200-residue output and exact longer consensus roundtrips
	 * through 1,000 residues; reject a 1,001-residue consensus. */
	private static void testConsensusLengthContract(Path root) throws Exception{
		final Fixture boundary=new Fixture(root.resolve("consensus_200"));
		boundary.consensus(200);
		final Path first=boundary.root.resolve("first"), second=boundary.root.resolve("second");
		FamilyReferenceExport.run(boundary.args(first));
		FamilyReferenceExport.run(boundary.args(second));
		check(Arrays.equals(Files.readAllBytes(first.resolve("references.tsv")), Files.readAllBytes(second.resolve("references.tsv"))),
			"Existing 200AA references remain byte-deterministic");
		for(int len:new int[]{241, 247, 1000}){
			final Fixture f=new Fixture(root.resolve("consensus_"+len));
			final String raw=sequence(len);
			f.consensus(raw);
			final Path out=f.root.resolve("out");
			FamilyReferenceExport.run(f.args(out));
			check(Arrays.equals(firstFamily(out).consensusRaw, raw.getBytes(StandardCharsets.US_ASCII)),
				"Long consensus raw bytes round-trip at "+len+"AA");
		}
		final Fixture tooLong=new Fixture(root.resolve("consensus_1001"));
		tooLong.consensus(1001);
		reject(tooLong.args(tooLong.root.resolve("out")));
	}
	/** Extend fresh fixture arguments with pinned complete alignment uppers and declared query count. */
	private static String[] alignmentArgs(Fixture f, Path out, Path table, int queries){
		final String[] base=f.args(out), args=Arrays.copyOf(base, base.length+3);
		args[base.length]="alignmentupper="+table;
		args[base.length+1]="alignmentuppersha80="+hash(table);
		args[base.length+2]="alignmentqueries="+queries;
		return args;
	}
	/** Load the first five-line record after a known five-line fixture header.
	 * This convenience reader assumes the exporter's output layout; it is not a file validator. */
	private static FamilyReference firstFamily(Path out) throws Exception{
		final List<String> lines=Files.readAllLines(out.resolve("references.tsv"), StandardCharsets.US_ASCII);
		final byte[][] record=new byte[5][];
		for(int i=0; i<5; i++){record[i]=lines.get(i+5).getBytes(StandardCharsets.US_ASCII);}
		return FamilyReferenceCodec.parseFiveLines(record, 0, new LineParser1('\t'));
	}

	/*--------------------------------------------------------------*/
	/*----------------     Fixture construction     ----------------*/
	/*--------------------------------------------------------------*/

	/** Owned mutable on-disk export inputs with one live and one empty family.
	 * Mutation cases explicitly rebind manifests when appropriate; generated files
	 * remain under the supplied test directory and are never automatically removed. */
	private static final class Fixture{
		/** Create a fresh fixture subtree, populate literal input tables and bind their manifests. */
		Fixture(Path root) throws Exception{
			this.root=root;
			Files.createDirectory(root);
			refs=Files.createDirectory(root.resolve("refs"));
			lower=Files.createDirectory(root.resolve("lower"));
			upper=Files.createDirectory(root.resolve("upper"));
			refIndex=root.resolve("references.tsv");
			upperIndex=root.resolve("uppers.tsv");
			consensus=root.resolve("consensus.faa");
			write(refs.resolve("members.tsv"),
				"roster_rank\tfamily_id\tmembers\tresidues\tvalid_words\n0\tf0\t1\t100\t1\n1\tempty\t0\t0\t0\n");
			final String row="tail\tf0\t1024\t1\t100\t1\t";
			write(refs.resolve("profiles.tsv"),
				DDLFamilyProfiles.HEADER16+"\n"+row+"literal\t"+array(1)+"\t"+ones()+"\t"+array(1)+"\n"+
				row+"nonempty\t"+array(1)+"\t"+array(1)+"\t"+array(1)+"\n");
			final ProteinSpectrum.Spectrum s=new ProteinSpectrum.Spectrum(2);
			s.counts[0]=s.total=99;
			final ByteBuilder spectrum=new ByteBuilder();
			ProteinSpectrum.writeRow(spectrum, "tail", "f0", s);
			write(refs.resolve("spectra.tsv"), ProteinSpectrum.header(2)+"\n"+spectrum);
			write(refs.resolve("summary.tsv"),
				"namespace\tk\tseed\texponent\tregisterbits\tbuckets\tshard\tbuilt_families\tzero_families\tmembers\tresidues\n"+
				"tail\t5\t0\t6\t16\t1024\t0\t1\t1\t1\t100\n");
			write(lower.resolve("cutoffs.tsv"),
				FamilyCalibration.CUTOFF_HEADER+"\n0\tf0\t1\tCALIBRATED\t0.2\t1\t30\t1\t-0.4\t1\n1\tempty\t0\tNO_MEMBERS\tNA\t0\tNA\t0\tNA\t0\n");
			write(lower.resolve("spectra.tsv"), "unused lower calibration spectra\n");
			write(upper.resolve("uppers.tsv"), FamilyGateExport.RAW_UPPER_HEADER+"\n0\tf0\t100000\t0.5\t500\t80\t-0.001\n");
			for(Path p:new Path[]{refs, lower, upper}){
				write(p.resolve("recipe.txt"), "fixture=true\n");
				write(p.resolve("timing.tsv"), "fixture\t0\n");
				if(!p.equals(refs)){write(p.resolve("summary.tsv"), "fixture\t0\n");}
			}
			writeUpperSummary(100000, 500, 500);
			write(refs.resolve("PASS"), "FAMILY_FULL_REFERENCES_PASS\n");
			write(lower.resolve("PASS"), "FAMILY_ALL_MEMBER_CALIBRATION_PASS\n");
			write(upper.resolve("PASS"), "FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS\n");
			consensus(100);
			bind();
		}
		/** Replace the consensus file with an all-alanine live sequence of the requested length. */
		void consensus(int len) throws Exception{
			final char[] aa=new char[len];
			Arrays.fill(aa, 'A');
			consensus(new String(aa));
		}
		/** Replace the live consensus bytes and retain the explicit empty-family placeholder. */
		void consensus(String sequence) throws Exception{write(consensus, ">f0\n"+sequence+"\n>empty\nX\n");}
		/** Refresh manifests and indices for the default one-member, 100-residue census. */
		void bind() throws Exception{
			bind(1, 100);
		}
		/** Refresh exact fixture manifests and index identities after deliberate input mutation.
		 * The caller supplies the expected live member/residue census; this does not recount data. */
		void bind(long members, long residues) throws Exception{
			manifest(refs, new String[]{"profiles.tsv", "spectra.tsv", "members.tsv", "summary.tsv", "recipe.txt",
				"timing.tsv", "PASS"});
			manifest(lower, new String[]{"cutoffs.tsv", "spectra.tsv", "summary.tsv", "recipe.txt", "timing.tsv", "PASS"});
			manifest(upper, new String[]{"uppers.tsv", "summary.tsv", "recipe.txt", "timing.tsv", "PASS"});
			final Path rm=refs.resolve("result.sha80.tsv"), lm=lower.resolve("result.sha80.tsv"),
				um=upper.resolve("result.sha80.tsv");
			write(refIndex, FamilyGateStats.INDEX_HEADERS[0]+"\n0\t2\t1\t1\t"+members+"\t"+residues+"\t"+rm+"\t"+hash(rm)+"\t"+lm+"\t"+hash(lm)+"\tNEW\n");
			write(upperIndex, FamilyGateStats.INDEX_HEADERS[1]+"\n0\t2\t1\t1\t"+members+"\t"+residues+"\t1\t"+um+"\t"+hash(um)+"\tNEW\n");
		}
		/** Write synthetic upper-summary counts with fixed timing fields for metadata checks. */
		void writeUpperSummary(int negatives, int top, int alignments) throws Exception{
			write(upper.resolve("summary.tsv"),
				"families\tprocessed_families\tnegatives\tcosine_calls\talignments\ttop\tthreads\tshard\tshards\tload_seconds\tddl_score_seconds\tcosine_seconds\talign_seconds\ttotal_seconds\n"+
				"1\t1\t"+negatives+"\t"+negatives+"\t"+alignments+"\t"+top+"\t1\t0\t1\t0\t0\t0\t0\t1\n");
		}
		/** Return fresh default-export arguments with current fixture index/consensus identities. */
		String[] args(Path out){
			return new String[]{"referenceindex="+refIndex, "referenceindexsha80="+hash(refIndex),
				"upperindex="+upperIndex, "upperindexsha80="+hash(upperIndex), "consensus="+consensus,
				"consensussha80="+hash(consensus), "out="+out, "maxfamilies=10"};
		}
		/** Owned root, component directories, index files and mutable consensus path. */
		final Path root, refs, lower, upper, refIndex, upperIndex, consensus;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Return 1,024 decimal registers: first followed by 1,023 zeroes. */
	private static String array(int first){
		final ByteBuilder b=new ByteBuilder().append(first);
		for(int i=1; i<1024; i++){b.append(",0");}
		return b.toString();
	}
	/** Return 1,024 comma-separated decimal ones for literal profile support. */
	private static String ones(){
		final ByteBuilder b=new ByteBuilder("1");
		for(int i=1; i<1024; i++){b.append(",1");}
		return b.toString();
	}
	/** Return deterministic ABZJX repeats to exercise raw consensus preservation. */
	private static String sequence(int len){
		final StringBuilder b=new StringBuilder(len);
		for(int i=0; i<len; i++){b.append("ABZJX".charAt(i%5));}
		return b.toString();
	}
	/** Replace the fixture manifest with identities for the explicitly named files in order. */
	private static void manifest(Path dir, String[] files) throws Exception{
		final ByteBuilder b=new ByteBuilder("file\tsha80\n");
		for(String f:files){b.append(f).tab().append(hash(dir.resolve(f))).nl();}
		write(dir.resolve("result.sha80.tsv"), b.toString());
	}
	/** Replace one zero-based TSV row/field and preserve every other field.
	 * Rewrites the synthetic ASCII file with LF endings and a terminal newline. */
	private static void cell(Path path, int row, int term, String value) throws Exception{
		final List<String> rows=Files.readAllLines(path, StandardCharsets.US_ASCII);
		final byte[] bytes=rows.get(row).getBytes(StandardCharsets.US_ASCII);
		final LineParser1 lp=new LineParser1('\t').set(bytes);
		lp.setBounds(term);
		rows.set(row, new String(bytes, 0, lp.a(), StandardCharsets.US_ASCII)+value+new String(bytes, lp.b(),
			bytes.length-lp.b(), StandardCharsets.US_ASCII));
		write(path, String.join("\n", rows)+"\n");
	}
	/** Read one zero-based ASCII fixture row, excluding its line terminator. */
	private static String line(Path p, int row) throws Exception{
		return Files.readAllLines(p, StandardCharsets.US_ASCII).get(row);
	}
	/** Replace the first matching key=value argument in place; fail if the fixture key is absent. */
	private static void replace(String[] args, String key, String value){
		for(int i=0; i<args.length; i++){
			if(args[i].startsWith(key+"=")){
				args[i]=key+"="+value;
				return;
			}
		}
		throw new AssertionError(key);
	}
	/** Require export to throw IllegalArgumentException and count that expected rejection. */
	private static void reject(String[] args) throws Exception{
		try{
			FamilyReferenceExport.run(args);
		}catch(IllegalArgumentException expected){
			checks++;
			return;
		}
		throw new AssertionError("Malformed export accepted");
	}
	/** Require the action to throw IllegalArgumentException; other failures propagate. */
	private static void reject(ThrowingAction action) throws Exception{
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			checks++;
			return;
		}
		throw new AssertionError("Malformed export accepted");
	}
	/** Require IllegalArgumentException containing the expected diagnostic fragment.
	 * A different exception or a normal return fails the case. */
	private static void rejectMessage(ThrowingAction action, String message) throws Exception{
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(message), "Expected rejection containing: "+message);
			return;
		}
		throw new AssertionError("Malformed export accepted");
	}
	/** Measure the fixture file identity required by legacy exporter input/output checks. */
	private static String hash(Path path){return DigestSuffix.file(path.toString());}
	/** Read a complete small uncompressed fixture as ASCII. */
	private static String read(Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	/** Read a complete small possibly compressed fixture as ASCII, closing the owned stream. */
	private static String readCompressed(Path path) throws Exception{
		final java.io.ByteArrayOutputStream out=new java.io.ByteArrayOutputStream();
		final java.io.InputStream in=ReadWrite.getInputStream(path.toString(), false, false);
		try{
			final byte[] buffer=new byte[4096];
			for(int len=in.read(buffer); len>=0; len=in.read(buffer)){out.write(buffer, 0, len);}
		}finally{in.close();}
		return new String(out.toByteArray(), StandardCharsets.US_ASCII);
	}
	/** Write or replace a synthetic ASCII fixture; the caller owns the destination. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	/** Fail immediately with the supplied reason, or increment the successful-check count. */
	private static void check(boolean yes, String message){
		if(!yes){throw new AssertionError(message);}
		checks++;
	}
	/** Fixture action that can propagate checked failures into an expected-rejection helper. */
	private interface ThrowingAction{
		/** Execute one fixture operation; unexpected checked failures propagate. */
		void run() throws Exception;
	}
	/*--------------------------------------------------------------*/
	/*----------------    Global state snapshot     ----------------*/
	/*--------------------------------------------------------------*/

	/** Snapshot of the Shared/ReadWrite settings touched by direct-BGZF test cases.
	 * Capture, mutate and restore serially; this is not a process-wide synchronization mechanism. */
	private static final class TestZipState{
		/** Store the captured scalar compression and shared-thread settings. */
		TestZipState(int threads_, int buffers_, boolean usePigz_, boolean useUnpigz_, boolean useBgzf_,
				boolean useBgzip_, boolean forceBgzip_, boolean allowNativeBgzf_, int zipLevel_, int zipThreads_,
				float zipThreadMult_, boolean setZipThreads_, boolean setZipThreadMult_){
			threads=threads_;
			buffers=buffers_;
			usePigz=usePigz_;
			useUnpigz=useUnpigz_;
			useBgzf=useBgzf_;
			useBgzip=useBgzip_;
			forceBgzip=forceBgzip_;
			allowNativeBgzf=allowNativeBgzf_;
			zipLevel=zipLevel_;
			zipThreads=zipThreads_;
			zipThreadMult=zipThreadMult_;
			setZipThreads=setZipThreads_;
			setZipThreadMult=setZipThreadMult_;
		}
		/** Snapshot the current effective values and explicit-setting flags. */
		static TestZipState capture(){
			return new TestZipState(Shared.threads(), Shared.READ_BUFFER_NUM_BUFFERS(),
				ReadWrite.USE_PIGZ, ReadWrite.USE_UNPIGZ, ReadWrite.USE_BGZF, ReadWrite.USE_BGZIP,
				ReadWrite.FORCE_BGZIP, ReadWrite.ALLOW_NATIVE_BGZF, ReadWrite.ZIPLEVEL,
				ReadWrite.MAX_ZIP_THREADS(), ReadWrite.ZIP_THREAD_MULT(), ReadWrite.SET_ZIP_THREADS,
				ReadWrite.SET_ZIP_THREAD_MULT);
		}
		/** Install deliberately nondefault globals so missing restoration becomes observable. */
		static void setNondefault(){
			Shared.setThreads(3);
			Shared.setBuffers(11);
			ReadWrite.USE_PIGZ=true;
			ReadWrite.USE_UNPIGZ=true;
			ReadWrite.USE_BGZF=false;
			ReadWrite.USE_BGZIP=false;
			ReadWrite.FORCE_BGZIP=false;
			ReadWrite.ALLOW_NATIVE_BGZF=true;
			ReadWrite.ZIPLEVEL=3;
			ReadWrite.setZipThreads(7);
			ReadWrite.setZipThreadMult(.5f);
			ReadWrite.SET_ZIP_THREADS=false;
			ReadWrite.SET_ZIP_THREAD_MULT=false;
		}
		/** Restore the snapshot, including explicit-setting flags after their setters run. */
		void restore(){
			Shared.setThreads(threads);
			Shared.setBuffers(buffers);
			ReadWrite.USE_PIGZ=usePigz;
			ReadWrite.USE_UNPIGZ=useUnpigz;
			ReadWrite.USE_BGZF=useBgzf;
			ReadWrite.USE_BGZIP=useBgzip;
			ReadWrite.FORCE_BGZIP=forceBgzip;
			ReadWrite.ALLOW_NATIVE_BGZF=allowNativeBgzf;
			ReadWrite.ZIPLEVEL=zipLevel;
			ReadWrite.setZipThreads(zipThreads);
			ReadWrite.setZipThreadMult(zipThreadMult);
			ReadWrite.SET_ZIP_THREADS=setZipThreads;
			ReadWrite.SET_ZIP_THREAD_MULT=setZipThreadMult;
		}
		/** Compare every captured setting with the current process globals. */
		boolean equalsCurrent(){
			return threads==Shared.threads() && buffers==Shared.READ_BUFFER_NUM_BUFFERS() &&
				usePigz==ReadWrite.USE_PIGZ && useUnpigz==ReadWrite.USE_UNPIGZ &&
				useBgzf==ReadWrite.USE_BGZF && useBgzip==ReadWrite.USE_BGZIP &&
				forceBgzip==ReadWrite.FORCE_BGZIP && allowNativeBgzf==ReadWrite.ALLOW_NATIVE_BGZF &&
				zipLevel==ReadWrite.ZIPLEVEL && zipThreads==ReadWrite.MAX_ZIP_THREADS() &&
				zipThreadMult==ReadWrite.ZIP_THREAD_MULT() &&
				setZipThreads==ReadWrite.SET_ZIP_THREADS &&
				setZipThreadMult==ReadWrite.SET_ZIP_THREAD_MULT;
		}
		/** Captured shared worker/buffer counts and compression level/thread limit. */
		final int threads, buffers, zipLevel, zipThreads;
		/** Captured compression thread multiplier. */
		final float zipThreadMult;
		/** Captured compressor/decompressor selection and native-BGZF flags. */
		final boolean usePigz, useUnpigz, useBgzf, useBgzip, forceBgzip, allowNativeBgzf;
		/** Whether the corresponding compression settings were explicitly selected. */
		final boolean setZipThreads, setZipThreadMult;
	}
	/** Process-local count of passed checks and expected rejections; not reset by main(). */
	private static int checks;
}
