package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import prot.DigestSuffix;
import prot.FamilyMemberScalars;
import prot.FamilyMemberScalarsCollect;
import structures.ByteBuilder;

/** Independent file-backed moments, sharding and worker-failure checks.
 * All fixture files remain under the caller's fresh Dori root.
 * @author Keqing
 */
public final class FamilyMemberScalarsTest{

	/** Compares two saved tables, or creates a fresh fixture root and checks export, collection and failure paths. */
	public static void main(String[] args) throws Exception{
		if(args.length==3 && args[0].startsWith("left=") && args[1].startsWith("right=") && args[2].startsWith("rows=")){
			final String[][] left=read(Paths.get(args[0].substring(5))), right=read(Paths.get(args[1].substring(6)));
			require(left.length==Integer.parseInt(args[2].substring(5)), "Wrong real-pilot family count");
			compare(left, right);
			System.err.println("FAMILY_MEMBER_SCALARS_PARITY_PASS rows="+left.length);
			return;
		}
		if(args.length!=1 || !args[0].startsWith("root=")){throw new IllegalArgumentException("Require root=fresh_fixture_directory");}
		final Path root=Paths.get(args[0].substring(5));
		Files.createDirectory(root);
		final Fixture fixture=new Fixture(root);
		final String[][] serial=run(fixture, "serial", 1, 1, 0);
		check(serial);
		final String[][] parallel=run(fixture, "parallel", 4, 1, 0);
		compare(serial, parallel);
		final String[][] even=run(fixture, "even", 2, 2, 0), odd=run(fixture, "odd", 2, 2, 1);
		require(even.length==2 && odd.length==2, "Shard roster sizes differ");
		compareRow(serial[0], even[0]);
		compareRow(serial[2], even[1]);
		compareRow(serial[1], odd[0]);
		compareRow(serial[3], odd[1]);
		hydrophobicity(fixture, serial);
		structural(fixture);
		collect(fixture, serial);
		// An admitted/pinned chunk with a wrong family identity must fail on a worker.
		final Path bad=fixture.files[0];
		write(bad, fixture.contents[0].replace("\tA\t", "\tWRONG\t"));
		fixture.manifest();
		fails(fixture, "bad_identity", "Member family rank/identity/bucket mismatch");
		write(bad, fixture.contents[0].replace(repeat('A', 100), repeat('X', 100)));
		fixture.manifest();
		fails(fixture, "bad_amino", "Accepted original tail contains a noncanonical amino acid");
		write(bad, fixture.contents[0]);
		fixture.manifest();
		moments();
		roundedStatistics();
		System.err.println("FAMILY_MEMBER_SCALARS_TEST_PASS checks="+checks);
	}
	/** Same literal member statistics without invented digests; metadata changes fail before publication. */
	private static void structural(final Fixture parent) throws Exception{
		final Path root=parent.root.resolve("structural").toAbsolutePath();
		Files.createDirectory(root);
		final Fixture f=new Fixture(root, true);
		final Path baseline=root.resolve("baseline.tsv");
		final String[] hydroSha80=hydroArguments(f, baseline, 1, 1, 0, true);
		final String[] sha80=java.util.Arrays.copyOf(hydroSha80, hydroSha80.length+1);
		sha80[sha80.length-1]="verification=sha80";
		FamilyMemberScalars.main(sha80);
		final StringBuilder manifest=new StringBuilder("batch_id\tbucket\tfile\tsha80\tmembers\tresidues\n");
		final List<String> rows=Files.readAllLines(f.chunks, StandardCharsets.US_ASCII);
		for(int i=1; i<rows.size(); i++){
			final String[] fields=rows.get(i).split("\t", -1);
			require(fields.length==6, "Fixture chunk schema changed");
			fields[3]="unmeasured";
			manifest.append(String.join("\t", fields)).append('\n');
		}
		write(f.chunks, manifest.toString());
		final Path out=root.resolve("result.tsv");
		FamilyMemberScalars.main(structuralArguments(f, out));
		require(java.util.Arrays.equals(Files.readAllBytes(baseline), Files.readAllBytes(out)),
			"Structural verification must preserve the same original-member KD bytes as SHA80 mode");
		final String[] legacy=structuralArguments(f, root.resolve("legacy_unknown.tsv"));
		legacy[legacy.length-1]="verification=sha80";
		rejectStructural(legacy, root.resolve("legacy_unknown.tsv"), "lowercase sha80");
		final String[] malformed=structuralArguments(f, root.resolve("bad_identity.tsv"));
		for(int i=0; i<malformed.length; i++){if(malformed[i].startsWith("familiessha80=")){malformed[i]="familiessha80=invalid";}}
		rejectStructural(malformed, root.resolve("bad_identity.tsv"), "Invalid structural member input identity");
		final String[] unknown=structuralArguments(f, root.resolve("bad_mode.tsv"));
		unknown[unknown.length-1]="verification=unknown";
		rejectStructural(unknown, root.resolve("bad_mode.tsv"), "verification=sha80|structural");
		final Path link=root.resolve("linked_chunk.tsv");
		Files.createSymbolicLink(link, f.files[0]);
		write(f.chunks, manifest.toString().replace(f.files[0].toString(), link.toString()));
		rejectStructural(structuralArguments(f, root.resolve("symlink.tsv")), root.resolve("symlink.tsv"), "absolute and nonsymlink");
		write(f.chunks, manifest.toString());
		final Path changed=root.resolve("changed.tsv");
		final FamilyMemberScalars prepared=new FamilyMemberScalars(structuralArguments(f, changed));
		write(f.files[0], f.contents[0]+"\n");
		// Invoke the existing private execution entry to mutate deterministically after input capture,
		// rather than relying on a race between a background writer and worker startup.
		final java.lang.reflect.Method process=FamilyMemberScalars.class.getDeclaredMethod("process");
		process.setAccessible(true);
		boolean rejected=false;
		try{process.invoke(prepared);}
		catch(java.lang.reflect.InvocationTargetException expected){rejected=hasReason(expected, "metadata changed");}
		require(rejected && !Files.exists(changed), "Changed structural input must fail before publishing metadata");
		write(f.files[0], f.contents[0]);
	}
	/** Builds the short-member, hydrophobicity-only structural command without measured input identities. */
	private static String[] structuralArguments(final Fixture f, final Path out){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>();
		args.add("families="+f.roster);
		args.add("familiessha80=unmeasured");
		args.add("chunks="+f.chunks);
		args.add("chunkssha80=unmeasured");
		args.add("out="+out);
		args.add("expectedfamilies=4");
		args.add("expectedmembers=6");
		args.add("batches=2");
		args.add("buckets=2");
		args.add("shards=1");
		args.add("shard=0");
		args.add("t=1");
		args.add("ow=f");
		args.add("scalars=hydrophobicity");
		args.add("minmemberaa=40");
		args.add("maxmemberaa=99");
		args.add("verification=structural");
		return args.toArray(new String[args.size()]);
	}
	/** Requires the structural exporter to explain rejection and leave its output absent. */
	private static void rejectStructural(final String[] args, final Path out, final String reason) throws Exception{
		boolean rejected=false;
		try{FamilyMemberScalars.main(args);}catch(Exception expected){rejected=hasReason(expected, reason);}
		require(rejected && !Files.exists(out), "Structural input rejection must explain the failure and publish nothing: "+reason);
	}
	/** Requires shard collection to explain rejection and leave its output absent. */
	private static void rejectCollect(final String[] args, final Path out, final String reason) throws Exception{
		boolean rejected=false;
		try{FamilyMemberScalarsCollect.main(args);}catch(Exception expected){rejected=hasReason(expected, reason);}
		require(rejected && !Files.exists(out), "Collect input rejection must explain the failure and publish nothing: "+reason);
	}
	/** Searches nested failure messages, including reflected and worker exceptions. */
	private static boolean hasReason(final Throwable failure, final String reason){
		for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
			if(cause.getMessage()!=null && cause.getMessage().contains(reason)){return true;}
		}
		return false;
	}
	/** Joins both saved shards under structural and exact verification and rejects missing exact identities. */
	private static void collect(final Fixture f, final String[][] expected)throws Exception{
		final Path even=f.root.resolve("even.tsv"), odd=f.root.resolve("odd.tsv");
		final Path structuralParts=f.root.resolve("parts_structural.tsv");
		final Path structuralOut=f.root.resolve("collected_structural.tsv");
		write(structuralParts, partsManifest(even, odd, false));
		FamilyMemberScalarsCollect.main(collectArguments(f, structuralParts, structuralOut,
			"unmeasured", "unmeasured", null));
		compare(expected, read(structuralOut));
		final Path exactParts=f.root.resolve("parts_exact.tsv"), exactOut=f.root.resolve("collected_exact.tsv");
		write(exactParts, partsManifest(even, odd, true));
		FamilyMemberScalarsCollect.main(collectArguments(f, exactParts, exactOut,
			DigestSuffix.file(f.roster.toString()), DigestSuffix.file(exactParts.toString()), "verification=sha80"));
		compare(expected, read(exactOut));
		final Path rejected=f.root.resolve("collect_rejected.tsv");
		rejectCollect(collectArguments(f, structuralParts, rejected, DigestSuffix.file(f.roster.toString()),
			DigestSuffix.file(structuralParts.toString()), "verification=sha80"), rejected, "lowercase sha80");
	}
	/** Describes the two ordered shard paths with either measured identities or explicit unmeasured markers. */
	private static String partsManifest(final Path even, final Path odd, final boolean measured){
		return "shard\tfile\tsha80\n0\t"+even+'\t'+(measured ? DigestSuffix.file(even.toString()) : "unmeasured")+
			"\n1\t"+odd+'\t'+(measured ? DigestSuffix.file(odd.toString()) : "unmeasured")+"\n";
	}
	/** Builds a four-family, six-member collection command with an optional verification override. */
	private static String[] collectArguments(final Fixture f, final Path parts, final Path out,
			final String familyPin, final String partsPin, final String extra){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>();
		args.add("families="+f.roster);
		args.add("familiessha80="+familyPin);
		args.add("parts="+parts);
		args.add("partssha80="+partsPin);
		args.add("out="+out);
		args.add("expectedfamilies=4");
		args.add("expectedmembers=6");
		args.add("shards=2");
		args.add("ow=f");
		if(extra!=null){args.add(extra);}
		return args.toArray(new String[args.size()]);
	}

	/** Literal member means, unequal peptide lengths, empty families, shards and independent worker reductions. */
	private static void hydrophobicity(final Fixture fixture, final String[][] all) throws Exception{
		final String[][] serial=runHydro(fixture, "hydro_serial", 1, 1, 0, false);
		checkHydro(serial);
		for(int i=0; i<all.length; i++){
			for(int j=0; j<3; j++){require(serial[i][j].equals(all[i][j]), "KD-only census must preserve default family identity/count");}
			if(i!=2){close(Double.parseDouble(all[i][3]), Double.parseDouble(serial[i][3]), 1e-12);}
		}
		compare(serial, runHydro(fixture, "hydro_parallel", 4, 1, 0, false));
		final Path shortRoot=fixture.root.resolve("short_hydro");
		Files.createDirectory(shortRoot);
		final Fixture shortFixture=new Fixture(shortRoot, true);
		final String[][] shortSerial=runHydro(shortFixture, "serial", 1, 1, 0, true);
		checkHydro(shortSerial);
		compare(shortSerial, runHydro(shortFixture, "parallel", 4, 1, 0, true));
		final double pooled=(1.8*40-3.5*80-3.9*99)/219;
		require(Math.abs(Double.parseDouble(shortSerial[0][3])-pooled)>.1,
			"Family KD must equally weight each member mean, not weight members by residue count");
		final String[][] even=runHydro(shortFixture, "even", 2, 2, 0, true), odd=runHydro(shortFixture, "odd", 2, 2, 1, true);
		compareRow(shortSerial[0], even[0]);
		compareRow(shortSerial[2], even[1]);
		compareRow(shortSerial[1], odd[0]);
		compareRow(shortSerial[3], odd[1]);
		write(shortFixture.files[0], shortFixture.contents[0].replace(repeat('A', 40), repeat('X', 40)));
		shortFixture.manifest();
		final Path bad=shortRoot.resolve("bad.tsv");
		boolean rejected=false;
		try{FamilyMemberScalars.main(hydroArguments(shortFixture, bad, 4, 1, 0, true));}
		catch(Exception e){
			for(Throwable cause=e; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains("noncanonical amino acid")){rejected=true;}
			}
		}
		require(rejected && !Files.exists(bad), "Invalid KD input must fail on its worker without publishing a table");
	}
	/** Exports one hydrophobicity shard and reads rows only after checking its exact four-column schema. */
	private static String[][] runHydro(Fixture f, String name, int threads, int shards, int shard, boolean shortRange)throws Exception{
		final Path out=f.root.resolve(name+".tsv");
		FamilyMemberScalars.main(hydroArguments(f, out, threads, shards, shard, shortRange));
		final List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
		require(lines.get(0).equals("family_rank\tfamily_id\tN\thydrophobicity_mean"), "KD-only export must declare the exact four-field schema");
		final String[][] rows=new String[lines.size()-1][];
		for(int i=1; i<lines.size(); i++){
			rows[i-1]=lines.get(i).split("\t", -1);
			require(rows[i-1].length==4, "KD-only rows cannot contain invented measurements of other scalars");
		}
		return rows;
	}
	/** Selects hydrophobicity alone, optionally restricting fixture members to 40–99 amino acids. */
	private static String[] hydroArguments(Fixture f, Path out, int threads, int shards, int shard, boolean shortRange){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>(java.util.Arrays.asList(arguments(f, out, threads, shards, shard)));
		args.add("scalars=hydrophobicity");
		if(shortRange){
			args.add("minmemberaa=40");
			args.add("maxmemberaa=99");
		}
		return args.toArray(new String[args.size()]);
	}
	/** Checks literal per-member means and the undefined value for the empty family. */
	private static void checkHydro(String[][] rows){
		require(rows.length==4, "KD-only output must retain all roster rows, including empty families");
		final long[] counts={3, 1, 0, 2};
		final double[] means={(1.8-3.5-3.9)/3, -.4, 0, .95};
		for(int i=0; i<rows.length; i++){
			require(rows[i][0].equals(""+i) && rows[i][1].equals(""+(char)('A'+i)) && rows[i][2].equals(""+counts[i]), "KD-only identity or member census differs");
			if(counts[i]==0){
				require(rows[i][3].equals("NA"), "An empty family must not receive an invented zero hydrophobicity");
			}else{close(means[i], Double.parseDouble(rows[i][3]), 1e-12);}
		}
	}

	/** Exports a complete scalar table or shard beneath the fixture root, then reads its checked-width rows. */
	private static String[][] run(final Fixture fixture, final String name, final int threads, final int shards, final int shard) throws Exception{
		final Path output=fixture.root.resolve(name+".tsv");
		FamilyMemberScalars.main(arguments(fixture, output, threads, shards, shard));
		return read(output);
	}
	/** Reads the 23-column table: three identity/count fields and four statistics for each of five scalars. */
	private static String[][] read(final Path output) throws Exception{
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		require(lines.get(0).split("\t").length==23, "Expected identity/count and twenty scalar fields");
		final String[][] result=new String[lines.size()-1][];
		for(int i=1; i<lines.size(); i++){
			result[i-1]=lines.get(i).split("\t", -1);
			require(result[i-1].length==23, "Wrong row width");
		}
		return result;
	}

	/** Builds the original six-member export command with freshly measured fixture identities and no overwrite. */
	private static String[] arguments(final Fixture f, final Path out, final int threads, final int shards, final int shard){
		return new String[]{"families="+f.roster, "familiessha80="+DigestSuffix.file(f.roster.toString()),
			"chunks="+f.chunks, "chunkssha80="+DigestSuffix.file(f.chunks.toString()), "out="+out,
			"expectedfamilies=4", "expectedmembers=6", "batches=2", "buckets=2", "shards="+shards,
			"shard="+shard, "t="+threads, "ow=f"};
	}

	/** Checks the synthetic roster and independently calculated scalar expectations. */
	private static void check(final String[][] rows){
		require(rows.length==4, "All four roster identities, including the empty family, must survive");
		final long[] n={3, 1, 0, 2};
		for(int i=0; i<4; i++){
			require(Integer.parseInt(rows[i][0])==i && rows[i][1].equals(""+(char)('A'+i)) &&
				Long.parseLong(rows[i][2])==n[i], "Family identity/count differs");
		}
		final double hydro=(1.8-3.5-3.9)/3;
		stat(rows[0], 0, hydro, Math.sqrt((1.8*1.8+3.5*3.5+3.9*3.9)/3-hydro*hydro), -3.9, 1.8);
		stat(rows[0], 1, 0, Math.sqrt(2.0/3), -1, 1);
		stat(rows[0], 2, 0, 0, 0, 0);
		stat(rows[0], 3, 1.0/3, Math.sqrt(2.0/9), 0, 1);
		stat(rows[1], 0, -0.4, 0, -0.4, -0.4);
		stat(rows[1], 1, 0, 0, 0, 0);
		stat(rows[1], 2, 0, 0, 0, 0);
		stat(rows[1], 3, 1, 0, 1, 1);
		stat(rows[1], 4, 6.015, 0, 6.015, 6.015);
		for(int i=3; i<23; i++){require(rows[2][i].equals("NA"), "Empty family statistic must be NA");}
		stat(rows[3], 0, 0.95, 1.85, -0.9, 2.8);
		stat(rows[3], 1, 0, 0, 0, 0);
		stat(rows[3], 2, 1, 0, 1, 1);
		stat(rows[3], 3, 0, 0, 0, 0);
		stat(rows[3], 4, 6.015, 0, 6.015, 6.015);
	}
	/** Checks one scalar's mean, population standard deviation, minimum and maximum to absolute tolerance. */
	private static void stat(String[] row, int scalar, double mean, double sd, double min, double max){
		final double[] expected={mean, sd, min, max};
		for(int j=0; j<4; j++){close(expected[j], Double.parseDouble(row[3+scalar*4+j]), 1e-9);}
	}
	/** Requires equal row counts and compares each corresponding family row. */
	private static void compare(final String[][] a, final String[][] b){
		require(a.length==b.length, "Parallel row count differs");
		for(int i=0; i<a.length; i++){compareRow(a[i], b[i]);}
	}
	/** Compares already width-checked rows: exact identities/NA markers, numeric values within 1e-9. */
	private static void compareRow(final String[] a, final String[] b){
		for(int i=0; i<a.length; i++){
			if(i<3 || a[i].equals("NA")){
				require(a[i].equals(b[i]), "Parallel/shard identity or undefined value differs");
			}else{close(Double.parseDouble(a[i]), Double.parseDouble(b[i]), 1e-9);}
		}
	}
	/** Requires a worker failure to reach main with its diagnostic and prevent output publication. */
	private static void fails(final Fixture f, final String label, final String diagnostic) throws Exception{
		final Path output=f.root.resolve(label+".tsv");
		boolean rejected=false;
		try{FamilyMemberScalars.main(arguments(f, output, 4, 1, 0));}
		catch(Exception e){
			for(Throwable cause=e; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(diagnostic)){rejected=true;}
			}
		}
		require(rejected && !Files.exists(output), "Worker failure must reach main and prevent output: "+label);
	}
	/** Checks merged moments near a large mean, where naive subtraction can lose variance. */
	private static void moments(){
		final FamilyMemberScalars.Moments a=new FamilyMemberScalars.Moments(1), b=new FamilyMemberScalars.Moments(1);
		final double[] values={1e9-1, 1e9-1, 1e9-1, 1e9-1, 1e9-1};
		a.add(0, values);
		java.util.Arrays.fill(values, 1e9);
		b.add(0, values);
		java.util.Arrays.fill(values, 1e9+1);
		b.add(0, values);
		a.merge(b);
		require(a.n[0]==3, "Merged count differs");
		for(int i=0; i<5; i++){
			close(1e9, a.mean[i], 1e-7);
			close(2, a.m2[i], 1e-7);
			close(1e9-1, a.min[i], 0);
			close(1e9+1, a.max[i], 0);
		}
	}
	/** Requires a finite observed value within an absolute tolerance of the independent expectation. */
	private static void close(double expected, double actual, double tolerance){
		require(Double.isFinite(actual) && Math.abs(actual-expected)<=tolerance, "Expected "+expected+" got "+actual);
	}
	/** Literal decimal expectations catch accidental use of unrounded moments or a different formatter. */
	private static void roundedStatistics(){
		final FamilyMemberScalars.Moments m=new FamilyMemberScalars.Moments(2);
		m.n[1]=30;
		for(int i=5; i<10; i++){
			m.mean[i]=0.12345678901234;
			m.m2[i]=30*0.11111111111114*0.11111111111114;
			m.min[i]=-0.22222222222224;
			m.max[i]=0.33333333333334;
		}
		final double[] out=new double[20];
		final ByteBuilder scratch=new ByteBuilder("discarded scratch");
		m.exportStatistics(1, out, scratch);
		final double[] expected={0.123456789012, 0.111111111111, -0.222222222222, 0.333333333333};
		for(int i=0; i<20; i++){require(Double.doubleToLongBits(out[i])==Double.doubleToLongBits(expected[i%4]), "Exact12-decimal statistic differs at "+i);}
		require(out[0]!=m.mean[5] && m.mean[5]==0.12345678901234, "Export rounds output without modifying accumulated moments");
		m.m2[5]=30e-28;
		m.exportStatistics(1, out, scratch);
		require(out[1]==0, "A positive standard deviation below table precision becomes exact zero");
		rejectStatistics(()->m.exportStatistics(0, out, scratch), "Empty family cannot fabricate measured zeros");
		rejectStatistics(()->m.exportStatistics(2, out, scratch), "Out-of-range family rejected");
		rejectStatistics(()->m.exportStatistics(1, new double[19], scratch), "Output needs all twenty fields");
		rejectStatistics(()->m.exportStatistics(1, out, null), "Formatting scratch is required");
		final FamilyMemberScalars.Moments aliased=new FamilyMemberScalars.Moments(4);
		aliased.n[0]=1;
		rejectStatistics(()->aliased.exportStatistics(0, aliased.mean, scratch), "Output cannot overwrite accumulated moments");
		m.m2[5]=-1;
		rejectStatistics(()->m.exportStatistics(1, out, scratch), "Negative variance rejected");
		m.m2[5]=0;
		m.mean[5]=Double.POSITIVE_INFINITY;
		rejectStatistics(()->m.exportStatistics(1, out, scratch), "Nonfinite moments cannot be formatted as measured statistics");
	}
	/** Requires argument validation to reject the proposed statistics export. */
	private static void rejectStatistics(Runnable action, String message){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError(message);
	}
	/** Counts an attempted check and fails immediately with its diagnostic. */
	private static void require(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}
	/** Attempted checks in this sequential test process; not reset between helper calls. */
	private static int checks;
	/** Replaces one synthetic ASCII fixture; original chunk contents are retained in Fixture. */
	private static void write(Path path, String text) throws Exception{Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	/** Builds a homopolymer of the requested residue and length. */
	private static String repeat(char c, int n){
		final char[] s=new char[n];
		java.util.Arrays.fill(s, c);
		return new String(s);
	}

	/** Four-family, six-member fixture owning its paths and the original text of four mutable chunks. */
	private static final class Fixture{
		/** Writes the default 100/200-residue members into an already-created root. */
		Fixture(Path root_) throws Exception{this(root_, false);}
		/** Writes the roster/chunks and their manifest, optionally using unequal 40–99-residue members. */
		Fixture(Path root_, boolean shortMembers) throws Exception{
			root=root_;
			roster=root.resolve("families.tsv");
			chunks=root.resolve("chunks.tsv");
			write(roster, "roster_rank\tfamily_id\tassigned\n0\tA\t3\n1\tB\t1\n2\tC\t0\n3\tD\t2\n");
			final String header="family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence\n";
			contents=new String[]{header+row(0, "A", "a", 0, 'A', shortMembers ? 40 : 100)+row(0, "A", "d", 1, 'D', shortMembers ? 80 : 100),
				header+row(1, "B", "g", 2, 'G', shortMembers ? 60 : 200)+row(3, "D", "f", 3, 'F', shortMembers ? 41 : 100),
				header+row(0, "A", "k", 4, 'K', shortMembers ? 99 : 100), header+row(3, "D", "w", 5, 'W', shortMembers ? 59 : 100)};
			for(int i=0; i<4; i++){
				files[i]=root.resolve("chunk"+i+".tsv");
				write(files[i], contents[i]);
			}
			manifest();
		}
		/** Rewrites counts and exact identities from current chunk bytes, including deliberate malformed cases. */
		void manifest() throws Exception{
			final StringBuilder s=new StringBuilder("batch_id\tbucket\tfile\tsha80\tmembers\tresidues\n");
			for(int i=0; i<4; i++){
				final List<String> rows=Files.readAllLines(files[i], StandardCharsets.US_ASCII);
				long residues=0;
				for(int j=1; j<rows.size(); j++){residues+=rows.get(j).split("\t", -1)[5].length();}
				s.append(i<2 ? "b00" : "b01").append('\t').append(i%2).append('\t').append(files[i]).append('\t')
					.append(DigestSuffix.file(files[i].toString())).append('\t').append(rows.size()-1).append('\t')
					.append(residues).append('\n');
			}
			write(chunks, s.toString());
		}
		/** Encodes one member with its family rank, synthetic genome/query identity and amino-acid sequence. */
		static String row(int rank, String family, String query, int ordinal, char aa, int len){
			return rank+"\t"+family+"\tg|"+query+"\tg\t"+ordinal+"\t"+repeat(aa, len)+"\n";
		}
		/** Caller-owned directory and the roster/manifest files written within it. */
		final Path root, roster, chunks;
		/** Two batches times two buckets, in manifest order; bytes may be changed by failure cases. */
		final Path[] files=new Path[4];
		/** Original chunk texts retained for deterministic restoration after mutation. */
		final String[] contents;
	}
}
