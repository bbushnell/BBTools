package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;

import prot.DigestSuffix;

/** Creates shared synthetic inputs for default-hash wrapper validation.
 * Paired parent construction has one public home in FamilyPairedEvalTest.
 * That setup enables process-wide amino-input mode; run in an isolated test JVM.
 * @author Brian, Ayaka, Fischl
 */
public final class FamilyDefaultHashInputFixture{

	/*--------------------------------------------------------------*/
	/*----------------        Input generation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates both input sets and freezes their manifests under one fresh absolute root.
	 * The parent directory must exist. This generates inputs and malformed variants;
	 * it does not execute the later default-hash wrapper comparisons.
	 */
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Require a fresh absolute fixture directory");}
		final Path supplied=Paths.get(args[0]);
		check(supplied.isAbsolute(), "Default-hash fixture root must be supplied as an absolute path");
		final Path root=supplied.toAbsolutePath().normalize();
		check(!Files.exists(root) && !Files.isSymbolicLink(root),
			"Default-hash fixture root must be fresh, absolute and nonsymlinked");
		Files.createDirectory(root);
		final Path pairedRoot=root.resolve("paired");
		Files.createDirectory(pairedRoot);
		final FamilyPairedEvalTest.ParentFixture pairedFixture=FamilyPairedEvalTest.prepareFixture(pairedRoot);
		final String[] constructionArgs=pairedFixture.constructionArgs;
		final FamilyPairedEvalTest.Inputs paired=FamilyPairedEvalTest.writeDefaultHashInputs(
			pairedRoot, pairedRoot.resolve("paired_eval"));
		final FamilyHitFloorEvalTest.Inputs floor=FamilyHitFloorEvalTest.writeDefaultHashInputs(root.resolve("floor"));
		final Path badUnsupported=FamilyHitFloorEvalTest.mutate(root.resolve("floor/bad_unsupported.tsv"),
			floor.floors, 6, 4, "5");
		final Path badInitial0=FamilyHitFloorEvalTest.mutate(root.resolve("floor/bad_initial0.tsv"),
			floor.floors, 6, 5, "1");
		final Path badWinnerBins=FamilyHitFloorEvalTest.mutate(root.resolve("floor/bad_winner_bins.tsv"),
			floor.floors, 6, 6, "1");
		final Path badDenominator=FamilyHitFloorEvalTest.mutate(root.resolve("floor/bad_denominator.tsv"),
			floor.floors, 4, 1, "0");
		final Path badRoster=FamilyHitFloorEvalTest.mutate(root.resolve("floor/bad_roster.tsv"),
			floor.floors, 6, 1, "other");
		final StringBuilder roles=new StringBuilder("#bbgap_default_hash_input_roles\t1\nrole\tpath\tsha80\n");
		final HashSet<Path> retained=new HashSet<Path>();
		file(roles, "PAIR_CONFIG", paired.config, root, retained);
		file(roles, "PAIR_MODEL", paired.model, root, retained);
		file(roles, "PAIR_MEDIANS", paired.medians, root, retained);
		file(roles, "PAIR_PANEL", paired.panel, root, retained);
		file(roles, "FLOOR_PARENT", floor.parent, root, retained);
		file(roles, "FLOOR_FITTING", floor.fitting, root, retained);
		file(roles, "FLOOR_TABLE", floor.floors, root, retained);
		file(roles, "FLOOR_TABLE_GZ", floor.compressed, root, retained);
		file(roles, "FLOOR_POS", floor.positive, root, retained);
		file(roles, "FLOOR_NEG", floor.negative, root, retained);
		file(roles, "FLOOR_BAD_UNSUPPORTED", badUnsupported, root, retained);
		file(roles, "FLOOR_BAD_INITIAL0", badInitial0, root, retained);
		file(roles, "FLOOR_BAD_WINNER_BINS", badWinnerBins, root, retained);
		file(roles, "FLOOR_BAD_DENOMINATOR", badDenominator, root, retained);
		file(roles, "FLOOR_BAD_ROSTER", badRoster, root, retained);
		final StringBuilder construction=new StringBuilder(
			"#bbgap_default_hash_construction_evidence\t1\nrole\tpath\tsha80\n");
		final Path packageRoot=pairedRoot.resolve("package");
		file(construction, "PAIR_PACKAGE_CONFIG", packageRoot.resolve("classifier.config"), root, retained);
		file(construction, "PAIR_PACKAGE_REFERENCES", packageRoot.resolve("references.tsv.gz"), root, retained);
		file(construction, "PAIR_PACKAGE_DDLPROFILES", packageRoot.resolve("ddlprofiles.tsv"), root, retained);
		file(construction, "PAIR_PACKAGE_PROFILE_PAYLOAD",
			packageRoot.resolve("profiles/profiles.tsv.gz"), root, retained);
		file(construction, "PAIR_PACKAGE_PROFILE_SUMMARY", packageRoot.resolve("profiles/summary.tsv"), root, retained);
		file(construction, "PAIR_RAW_CALIBRATION", packageRoot.resolve("raw_calibration.tsv.gz"), root, retained);
		file(construction, "PAIR_PREPARATION_TIMINGS", packageRoot.resolve("preparation.tsv"), root, retained);
		file(construction, "PAIR_PREPARATION_INPUTS", packageRoot.resolve("preparation_inputs.tsv"), root, retained);
		file(construction, "PAIR_SOURCE_LENGTHS", pairedFixture.lengths, root, retained);
		file(construction, "PAIR_CALIBRATION_SUMMARY", packageRoot.resolve("summary.tsv"), root, retained);
		file(construction, "PAIR_CALIBRATION_FAMILIES",
			packageRoot.resolve("family_false_accepts.tsv.gz"), root, retained);
		file(construction, "PAIR_CALIBRATION_ROUNDS", packageRoot.resolve("rounds.tsv"), root, retained);
		file(construction, "PAIR_CALIBRATION_COHORTS", packageRoot.resolve("cohorts.tsv"), root, retained);
		file(construction, "PAIR_CALIBRATION_TABLE", packageRoot.resolve("calibration.tsv"), root, retained);
		file(construction, "PAIR_CALIBRATION_INPUTS", packageRoot.resolve("inputs.tsv"), root, retained);
		file(construction, "PAIR_CALIBRATION_INTERPRETATION",
			packageRoot.resolve("interpretation.txt"), root, retained);
		optionalFile(construction, "PAIR_CALIBRATION_FROZEN", packageRoot.resolve("FROZEN"), root, retained);
		file(construction, "PAIR_CALIBRATION_PASS", packageRoot.resolve("PASS"), root, retained);
		constructionInputs(construction, packageRoot.resolve("preparation_inputs.tsv"), root, retained);
		final StringBuilder dependencies=new StringBuilder(
			"#bbgap_default_hash_input_dependencies\t1\nsource\tpath\tsha80\n");
		final HashSet<Path> seen=new HashSet<Path>();
		dependencies(dependencies, seen, retained, root, "pair", paired.config);
		dependencies(dependencies, seen, retained, root, "floor", floor.parent);
		final Path manifest=root.resolve("manifest.tsv"), dependencyTable=root.resolve("dependencies.tsv");
		final Path generationTable=root.resolve("generation.tsv"), constructionTable=root.resolve("construction.tsv");
		Files.write(manifest, roles.toString().getBytes(StandardCharsets.US_ASCII),
			java.nio.file.StandardOpenOption.CREATE_NEW);
		Files.write(constructionTable, construction.toString().getBytes(StandardCharsets.US_ASCII),
			java.nio.file.StandardOpenOption.CREATE_NEW);
		Files.write(dependencyTable, dependencies.toString().getBytes(StandardCharsets.US_ASCII),
			java.nio.file.StandardOpenOption.CREATE_NEW);
		Files.write(generationTable, generation(constructionArgs).getBytes(StandardCharsets.US_ASCII),
			java.nio.file.StandardOpenOption.CREATE_NEW);
		retain(root, retained, manifest);
		retain(root, retained, constructionTable);
		retain(root, retained, dependencyTable);
		retain(root, retained, generationTable);
		makeReadonly(root, retained);
		System.out.println("FAMILY_DEFAULT_HASH_INPUT_FIXTURE_PASS manifest="+root.resolve("manifest.tsv"));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Manifest records      ----------------*/
	/*--------------------------------------------------------------*/

	/** Records the actual public producer and the arguments it passed to calibration. */
	private static String generation(String[] constructionArgs){
		final StringBuilder sb=new StringBuilder(
			"#bbgap_default_hash_input_generation\t1\nrecord\tphase\tordinal\tvalue\n");
		sb.append("tool\tpaired_parent\t0\tFamilyCalibrate\n");
		sb.append("source\tpaired_parent\t0\tFamilyPairedEvalTest.ParentFixture\n");
		for(int i=0; i<constructionArgs.length; i++){
			sb.append("argument\tpaired_parent\t").append(i).append('\t').append(constructionArgs[i]).append('\n');
		}
		return sb.toString();
	}
	/** Records supported config resources and their profile payload/summary dependencies.
	 * Relative paths resolve against their declaring file; seen paths are emitted
	 * once across both source packages, while every path remains in retained.
	 */
	private static void dependencies(StringBuilder sb, HashSet<Path> seen, HashSet<Path> retained, Path root,
			String source, Path config)throws Exception{
		file(sb, seen, retained, root, source, config);
		for(String line:Files.readAllLines(config, StandardCharsets.US_ASCII)){
			if(line.isEmpty() || line.charAt(0)=='#' || line.indexOf('=')<0){continue;}
			final String[] pair=line.split("=", 2);
			if(pair[0].equals("familyreferences") || pair[0].equals("ddlprofiles") || pair[0].equals("scalarbounds") ||
					pair[0].equals("blacklist") || pair[0].equals("lengthmedians") || pair[0].equals("familyhydro")){
				final Path declared=Paths.get(pair[1]);
				final Path path=(declared.isAbsolute() ? declared : config.getParent().resolve(declared)).
					toAbsolutePath().normalize();
				file(sb, seen, retained, root, source, path);
				if(pair[0].equals("ddlprofiles")){
					final Path base=path.getParent();
					for(String row:Files.readAllLines(path, StandardCharsets.US_ASCII)){
						if(row.trim().isEmpty() || row.charAt(0)=='#'){continue;}
						final String profile=row.split("\t", -1)[0];
						final Path payload=base.resolve(profile).toAbsolutePath().normalize();
						file(sb, seen, retained, root, source, payload);
						file(sb, seen, retained, root, source,
							payload.getParent().resolve("summary.tsv").toAbsolutePath().normalize());
					}
				}
			}
		}
	}
	/** Copies declared preparation inputs into the construction manifest after validating each row. */
	private static void constructionInputs(StringBuilder sb, Path table, Path root,
		HashSet<Path> retained)throws Exception{
		for(String line:Files.readAllLines(table, StandardCharsets.US_ASCII)){
			if(line.equals("path\tsha80")){continue;}
			final String[] fields=line.split("\t", -1);
			check(fields.length==2, "Malformed preparation input row: "+line);
			file(sb, "PAIR_PREPARATION_INPUT", Paths.get(fields[0]), root, retained);
		}
	}
	/** Records an optional artifact only when it exists; normal file validation still applies. */
	private static void optionalFile(StringBuilder sb, String role, Path path, Path root,
			HashSet<Path> retained)throws Exception{
		if(Files.exists(path, LinkOption.NOFOLLOW_LINKS)){file(sb, role, path, root, retained);}
	}
	/** Appends a role/path/digest row and retains its checked file for final read-only handling. */
	private static void file(StringBuilder sb, String role, Path path, Path root,
		HashSet<Path> retained)throws Exception{
		path=retain(root, retained, path);
		sb.append(role).append('\t').append(path).append('\t').append(DigestSuffix.file(path.toString())).append('\n');
	}
	/** Retains a checked dependency and emits its source/path/digest row only on first occurrence. */
	private static void file(StringBuilder sb, HashSet<Path> seen, HashSet<Path> retained, Path root,
			String source, Path path)throws Exception{
		path=retain(root, retained, path);
		if(seen.add(path)){
			sb.append(source).append('\t').append(path).append('\t')
				.append(DigestSuffix.file(path.toString())).append('\n');
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Path ownership        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a normalized regular file whose real path stays within the fixture root.
	 * The file itself must not be a symlink. Returns the normalized path used
	 * by the manifests; parent components are resolved only for containment checks.
	 */
	private static Path retain(Path root, HashSet<Path> retained, Path path)throws Exception{
		path=path.toAbsolutePath().normalize();
		final Path realRoot=root.toRealPath();
		check(Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS) && !Files.isSymbolicLink(path),
			"Retained fixture path must be a regular nonsymlink file: "+path);
		check(path.toRealPath().startsWith(realRoot),
			"Retained fixture path escaped the fixture root: "+path);
		retained.add(path);
		return path;
	}
	/** Rechecks every retained file before changing permissions, then removes write access.
	 * A failed permission change aborts; earlier successful changes are not rolled back.
	 */
	private static void makeReadonly(Path root, HashSet<Path> retained)throws Exception{
		final Path realRoot=root.toRealPath();
		for(Path path:retained){
			check(Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS) && !Files.isSymbolicLink(path),
				"Final fixture path must be a regular nonsymlink file: "+path);
			check(path.toRealPath().startsWith(realRoot),
				"Final fixture path escaped the fixture root: "+path);
		}
		for(Path path:retained){check(path.toFile().setReadOnly(), "Could not finalize readonly file: "+path);}
	}
	/** Rejects malformed fixture inputs or violated file-ownership contracts. */
	private static void check(boolean valid, String why){if(!valid){throw new IllegalArgumentException(why);}}
}
