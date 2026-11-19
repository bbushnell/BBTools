package prot;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.bam.BgzfSettings;

/** Collects disjoint FamilyMemberScalars shards without recomputing moments.
 * Keeps only aggregate rows (not members) in a dense rank-indexed array. Checks
 * every input before opening output, then emits original row bytes plus LF.
 * Single-reader ByteFile/ByteStreamWriter pattern: parallel parsing is needless
 * for this small aggregate table. Exact mode pins compressed parts as stored bytes.
 * @author Sayu
 */
public final class FamilyMemberScalarsCollect{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilyMemberScalarsCollect.class, false);
		try{new FamilyMemberScalarsCollect(pp.args).process();}finally{Shared.closeStream(pp.outstream);}
	}

	private FamilyMemberScalarsCollect(final String[] args) throws Exception{
		assert(args!=null) : "Collector needs explicit pinned input arguments, not implicit discovery";
		final Parser parser=new Parser();
		for(final String arg : args){
			final String[] kv=Parse.splitOnFirst(arg, '=');
			final String k=kv[0].toLowerCase(Locale.ROOT), v=kv.length>1 ? kv[1] : null;
			if(k.equals("families")){families=v;}
			else if(k.equals("familiessha80")){familyPin=v;}
			else if(k.equals("parts")){manifest=v;}
			else if(k.equals("partssha80")){partsPin=v;}
			else if(k.equals("expectedfamilies")){nf=Integer.parseInt(v);}
			else if(k.equals("expectedmembers")){expectedMembers=Long.parseLong(v);}
			else if(k.equals("shards")){shards=Integer.parseInt(v);}
			else if(k.equals("verification")){
				require(!verificationSeen && ("sha80".equals(v) || "structural".equals(v)),
					"Specify verification=sha80|structural once");
				verificationSeen=true; structural="structural".equals(v);
			}
			else if(!parser.parse(arg, k, v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		require(families!=null && manifest!=null && parser.out1!=null && !parser.append && !parser.overwrite,
			"Require families= parts= out= ow=f; append and overwrite are forbidden");
		require(nf>0 && expectedMembers>=0 && shards>0 && shards<=nf, "Invalid collection dimensions");
		out=parser.out1;
		require(new File(out).isAbsolute(), "Output must be an absolute regular-file path");
		destination=new File(out).getCanonicalFile();
		fresh();
		Shared.setThreads(1);
		// Aggregate-only parsing: avoid the fast decimal parser's long overflow
		// on malformed oversized fields. Member payloads never pass through here.
		Tools.FORCE_JAVA_PARSE_DOUBLE=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		BgzfSettings.USE_MULTITHREADED_BGZF=false;
		ids=new String[nf]; expected=new long[nf]; rows=new byte[nf][];
	}

	private void process() throws Exception{
		assert(rows.length==nf && expected.length==nf) : "Dense output indices must share roster dimensions";
		readRoster(); readManifest();
		int count=0; long members=0;
		for(final Part part : parts){readPart(part);}
		for(int rank=0; rank<nf; rank++){
			require(rows[rank]!=null, "Missing family rank: "+rank);
			count++; members=Math.addExact(members, expected[rank]);
		}
		require(count==nf && members==expectedMembers, "Collected population differs from required totals");
		verifyInputs();
		fresh();
		final ByteStreamWriter writer=new ByteStreamWriter(out, false, false, false);
		writer.start();
		try{
			writer.println(HEADER);
			for(final byte[] row : rows){writer.println(row);}
		}finally{require(!writer.poisonAndWait(), "Collection output failed");}
		verifyInputs();
		System.err.println("FAMILY_MEMBER_SCALARS_COLLECT_PASS families="+count+" members="+members+
			" shards="+shards+" original_rows=true");
	}

	private void readRoster() throws Exception{
		assert(ids.length==nf) : "Each dense roster rank must have one identity slot";
		rosterState=verifyInput(families, familyPin);
		final ByteFile input=ByteFile.makeByteFile1(families, false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> seen=new HashSet<String>();
		int rank=0; long members=0;
		try{
			require(Arrays.equals(input.nextLine(), ascii("roster_rank\tfamily_id\tassigned")), "Wrong family roster header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==3 && rank<nf && integer(lp, line, 0)==rank, "Roster must be dense and ordered");
				final String id=lp.parseString(1); final long n=integer(lp, line, 2);
				require(!id.isEmpty() && seen.add(id), "Duplicate or empty roster family identity");
				ids[rank]=id; expected[rank]=n; members=Math.addExact(members, n); rank++;
			}
		}finally{require(!input.close(), "Roster read failed");}
		require(rank==nf && members==expectedMembers, "Roster population differs from required totals");
	}

	/** Exactly one declared identity and canonical path for each shard. */
	private void readManifest() throws Exception{
		assert(parts.isEmpty()) : "A collector invocation must not accumulate manifests from prior runs";
		manifestState=verifyInput(manifest, partsPin);
		final ByteFile input=ByteFile.makeByteFile1(manifest, false);
		final LineParser1 lp=new LineParser1('\t');
		final boolean[] seen=new boolean[shards];
		final HashSet<File> paths=new HashSet<File>();
		try{
			require(Arrays.equals(input.nextLine(), ascii("shard\tfile\tsha80")), "Wrong parts manifest header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==3, "Part manifest must have three columns");
				final long number=integer(lp, line, 0);
				require(number<shards && !seen[(int)number], "Duplicate or out-of-range part shard");
				final String path=lp.parseString(1), hash=lp.parseString(2);
				final File file=inputPath(path);
				require(paths.add(file) && !file.equals(new File(families).getCanonicalFile()) &&
					!file.equals(new File(manifest).getCanonicalFile()), "Duplicate or control-alias part path");
				require(hash.matches("[0-9a-f]{20}") || structural && hash.equals("unmeasured"),
					"Part identity must be a lowercase sha80, or unmeasured in structural mode");
				parts.add(new Part((int)number, path, hash, structural ? new FileState(new File(path).toPath()) : null));
				seen[(int)number]=true;
			}
		}finally{require(!input.close(), "Part manifest read failed");}
		require(parts.size()==shards, "Missing part shard");
	}

	private void readPart(final Part part){
		assert(part.shard>=0 && part.shard<shards) : "Manifest validation supplies the rank modulo ownership check";
		if(part.state==null){pin(part.path, part.hash);}else{part.state.verify();}
		final ByteFile input=ByteFile.makeByteFile1(part.path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(input.nextLine(), HEADER), "Wrong 23-column member scalar header: "+part.path);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==23, "Member scalar row must have 23 columns");
				final long number=integer(lp, line, 0);
				require(number<nf, "Family rank is outside roster");
				final int rank=(int)number;
				require(rank%shards==part.shard, "Family rank is in the wrong shard: "+rank);
				require(rows[rank]==null, "Duplicate family rank: "+rank);
				require(lp.termEquals(ids[rank], 1), "Family identity differs at rank: "+rank);
				final long n=integer(lp, line, 2);
				require(n==expected[rank], "Family N differs at rank: "+rank);
				for(int scalar=0; scalar<5; scalar++){
					final int column=3+4*scalar;
					if(n==0){
						for(int j=0; j<4; j++){require(lp.termEquals("NA", column+j), "N=0 statistics must all be NA");}
					}else{
						final double mean=decimal(lp, line, column), sd=decimal(lp, line, column+1);
						final double min=decimal(lp, line, column+2), max=decimal(lp, line, column+3);
						require(sd>=0 && min<=mean && mean<=max, "Invalid scalar SD/range at rank: "+rank);
						require(n!=1 || (sd==0 && min==mean && mean==max), "N=1 requires SD=0 and equal extrema");
					}
				}
				rows[rank]=line;
			}
		}finally{require(!input.close(), "Part read failed: "+part.path);}
		if(part.state!=null){part.state.verify();}
	}

	/** Strict unsigned decimal count parsing; detects overflow before accepting it. */
	private static long integer(final LineParser1 lp, final byte[] line, final int column){
		assert(column>=0 && column<lp.terms()) : "Schema width must be checked before parsing a count";
		lp.setBounds(column); require(lp.b()>lp.a(), "Empty integer field");
		long value=0;
		for(int i=lp.a(); i<lp.b(); i++){
			final int digit=line[i]-'0'; require(digit>=0 && digit<=9, "Malformed unsigned integer field");
			value=Math.addExact(Math.multiplyExact(value, 10), digit);
		}
		return value;
	}

	/** Producer emits fixed decimals, never exponents, suffixes or whitespace. */
	private static double decimal(final LineParser1 lp, final byte[] line, final int column){
		assert(column>=3 && column<23) : "Only scalar fields use the producer's fixed-decimal format";
		lp.setBounds(column); int i=lp.a(), digits=0, dots=0;
		if(i<lp.b() && line[i]=='-'){i++;}
		for(; i<lp.b(); i++){
			if(line[i]=='.'){dots++;}
			else{require(line[i]>='0' && line[i]<='9', "Malformed/nonfinite scalar field"); digits++;}
		}
		require(digits>0 && dots<=1, "Malformed scalar decimal");
		final double value=lp.parseDouble(column);
		require(Double.isFinite(value), "Nonfinite scalar field"); return value;
	}

	private File inputPath(final String path) throws Exception{
		assert(destination!=null) : "Input alias checks require the canonical destination";
		final File file=new File(path);
		require(file.isAbsolute() && file.isFile(), "Input must be an absolute existing regular file: "+path);
		final File canonical=file.getCanonicalFile();
		require(!canonical.equals(destination), "Output aliases an input: "+path); return canonical;
	}
	private void fresh(){
		require(!Files.exists(new File(out).toPath(), LinkOption.NOFOLLOW_LINKS), "Output exists; use a fresh result path");
	}
	private FileState verifyInput(final String path, final String identity) throws Exception{
		inputPath(path);
		if(!structural){pin(path, identity); return null;}
		require(identity!=null && (identity.matches("[0-9a-f]{20}") || identity.equals("unmeasured")),
			"Invalid structural member input identity: "+path);
		return new FileState(new File(path).toPath());
	}
	/** Recheck all declared inputs before and after publication without hashing in structural mode. */
	private void verifyInputs(){
		if(!structural){
			pin(families, familyPin);
			pin(manifest, partsPin);
			for(final Part part : parts){pin(part.path, part.hash);}
			return;
		}
		rosterState.verify();
		manifestState.verify();
		for(final Part part : parts){part.state.verify();}
	}
	/** Immutable metadata captured before reading a file; never a content hash. */
	private static final class FileState{
		FileState(final Path path_){
			path=path_;
			initial=read();
		}
		private BasicFileAttributes read(){
			require(path.isAbsolute() && !Files.isSymbolicLink(path),
				"Structural member input must be absolute and nonsymlink: "+path);
			try{
				final BasicFileAttributes state=
					Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
				require(state.isRegularFile(), "Structural member input must be a regular file: "+path);
				return state;
			}catch(java.io.IOException e){
				throw new IllegalArgumentException("Cannot inspect structural member input: "+path, e);
			}
		}
		void verify(){
			final BasicFileAttributes current=read();
			require(initial.size()==current.size() &&
					initial.lastModifiedTime().equals(current.lastModifiedTime()) &&
					java.util.Objects.equals(initial.fileKey(), current.fileKey()),
				"Structural member input metadata changed: "+path);
		}
		final Path path;
		final BasicFileAttributes initial;
	}
	private static void pin(final String path, final String hash){
		require(hash!=null && hash.matches("[0-9a-f]{20}"), "Input pin must be a lowercase sha80");
		require(DigestSuffix.file(path).equals(hash), "Input content differs: "+path);
	}
	private static void require(final boolean ok, final String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static byte[] ascii(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static final class Part{
		Part(final int s, final String p, final String h, final FileState state_){
			shard=s; path=p; hash=h; state=state_;
		}
		final int shard; final String path, hash;
		final FileState state;
	}
	private static final byte[] HEADER=ascii("family_rank\tfamily_id\tN"+
		"\thydrophobicity_mean\thydrophobicity_sd\thydrophobicity_min\thydrophobicity_max"+
		"\tnet_charge_mean\tnet_charge_sd\tnet_charge_min\tnet_charge_max"+
		"\taromaticity_mean\taromaticity_sd\taromaticity_min\taromaticity_max"+
		"\tfraction_small_mean\tfraction_small_sd\tfraction_small_min\tfraction_small_max"+
		"\tpI_tail_mean\tpI_tail_sd\tpI_tail_min\tpI_tail_max");
	private String families, familyPin, manifest, partsPin, out;
	private boolean structural=true, verificationSeen;
	private FileState rosterState, manifestState;
	private int nf=79565, shards=64;
	private long expectedMembers=124056489;
	private final File destination;
	private final String[] ids;
	private final long[] expected;
	private final byte[][] rows;
	private final ArrayList<Part> parts=new ArrayList<Part>();
}
