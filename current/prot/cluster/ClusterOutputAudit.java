package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import parse.LineParser1;
import prot.DigestSuffix;
import prot.FamilyCalibration.VerificationMode;
import structures.IntHashMap;

/** Read-only conservation audit of an uncalibrated cluster export.
 * Streams assignments with O(families) counters; the canonical codec validates
 * model payloads. Does not score sequences, verify FASTA IDs, or establish
 * biological correctness. SHA80 mode checks pinned stored bytes before/after
 * inspection. Structural mode checks existence and parsed conservation only;
 * supplied identities are unmeasured claims, not content-integrity evidence.
 * @author Sayu, Fischl
 */
public final class ClusterOutputAudit{

	/** Static audit utility; no instance state. */
	private ClusterOutputAudit(){}

	/*--------------------------------------------------------------*/
	/*----------------         Entry points         ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute the complete audit before writing its receipt to stdout.
	 * Input-check failures print no success receipt; stdout failures may occur
	 * after partial output and are reported separately. */
	public static void main(final String[] args){
		final String report=audit(args);
		System.out.print(report);
		if(System.out.checkError()){throw new IllegalStateException("Could not write cluster audit receipt");}
	}

	/** Returns a complete receipt only after every check; never creates or changes input files. */
	public static String audit(final String[] args){
		final HashMap<String, String> o=options(args);
		final VerificationMode verification=verification(o);
		final HashMap<String, String> identities=new HashMap<String, String>();
		for(String key:FILES){identities.put(key, verify(o, key, verification));}
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(required(o, "codec"));
		final ClusterSketchCodec.Header h=loaded.header;
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(integer(o, "k"), required(o, "aabits"));
		final ProteinPositionMode position=ProteinPositionMode.parse(required(o, "positionmode"));
		final String blacklist=required(o, "blacklistsha80");
		require(blacklist.equals("none") || blacklist.matches("[0-9a-f]{20}") ||
			verification==VerificationMode.STRUCTURAL && blacklist.equals("unmeasured"),
			"Invalid blacklist identity for the verification policy");
		require(h.geometry.equals(geometry) && h.k==integer(o, "k"), "Codec k/alphabet differs from required geometry");
		require(h.width==integer(o, "width") && h.sketchSeed==Long.parseLong(required(o, "sketchseed")),
			"Codec width/sketch seed differs");
		require(h.positionMode==position && h.partitions==integer(o, "partitions") && h.overlap==integer(o, "overlap"),
			"Codec positional layout differs");
		require(Objects.equals(h.blacklistSha80, blacklist.equals("none") ? null : blacklist), "Codec blacklist pin differs");
		final HashMap<String, String> stats=stats(required(o, "stats"));
		final long expected=number(o, "inputs");
		require(number(stats, "input_records")==expected, "Stats input count differs from declared population");
		require(number(stats, "families")==loaded.families.length, "Stats family count differs from codec");
		require(ProteinKmerGeometry.make(integer(stats, "k"), required(stats, "aabits")).equals(geometry),
			"Stats k/alphabet differs");
		for(String key:new String[]{"width", "partitions", "overlap", "minlen", "maxalignlen"}){
			require(number(stats, key)==number(o, key), "Stats differs for "+key);
		}
		require(Long.parseLong(required(stats, "sketchseed"))==h.sketchSeed, "Stats sketch seed differs");
		require(required(stats, "positionmode").equals(position.flag) &&
			required(stats, "positional").equals(Boolean.toString(position.positional())), "Stats positional layout differs");
		require(required(stats, "blacklist_sha80").equals(blacklist), "Stats blacklist differs");
		final String ranking=required(o, "memberranking");
		require(ranking.equals("true") || ranking.equals("false"), "memberranking must be true or false");
		require(required(stats, "memberranking").equals(ranking), "Stats member ranking differs");
		require(number(stats, "completed_passes")<=number(stats, "passes"), "Completed passes exceed requested passes");
		final long[] counts=members(required(o, "members"), loaded, expected);
		require(number(stats, "assigned_records")==counts[0], "Stats assigned count differs from member rows");
		require(number(stats, "unassigned_records")==counts[1], "Stats unassigned count differs from member rows");
		require(number(stats, "filtered_records")==counts[2], "Stats filtered count differs from member rows");
		require(number(stats, "eligible_records")==counts[0]+counts[1],
			"Eligible must equal assigned plus unassigned (ClusteringOutput.summary)");
		assert(counts[0]+counts[1]+counts[2]==expected) : "Each member row has exactly one status; protects conservation receipt";
		for(String key:FILES){verify(o, key, verification);}
		final StringBuilder report=new StringBuilder("key\tvalue\n");
		for(String key:FILES){report.append(key).append("_sha80\t").append(identities.get(key)).append('\n');}
		report.append("geometry\t").append(geometry.name()).append("\nwidth\t").append(h.width)
			.append("\ninput_records\t").append(expected).append("\nfamilies\t").append(loaded.families.length)
			.append("\nassigned_records\t").append(counts[0]).append("\nunassigned_records\t").append(counts[1])
			.append("\nfiltered_records\t").append(counts[2]).append("\nstatus\tCLUSTER_OUTPUT_AUDIT_PASS\n");
		return report.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------   Assignment conservation    ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams every input ordinal once, checking each stable family against its codec count. */
	private static long[] members(final String path, final ClusterSketchCodec.Loaded loaded, final long expected){
		assert(loaded.header.families==loaded.families.length) : "Codec.read validates declared roster; counter ranks must agree";
		final IntHashMap ranks=new IntHashMap(Math.max(16, loaded.families.length));
		final long[] observed=new long[loaded.families.length];
		for(int i=0; i<observed.length; i++){ranks.put(loaded.families[i].familyId, i);}
		final ByteFile input=ByteFile.makeByteFile1(path, false);
		final LineParser1 lp=new LineParser1('\t');
		long ordinal=0, assigned=0, unassigned=0, filtered=0;
		try{
			require(Arrays.equals(input.nextLine(), MEMBER_HEADER), "Wrong clustering member header");
			for(byte[] line; (line=input.nextLine())!=null; ordinal++){
				lp.set(line);
				require(lp.terms()==4, "Member table requires four columns");
				require(ordinal<expected && unsignedField(lp, 0)==ordinal,
					"Member ordinal missing, duplicated, out of order, or beyond population");
				require(!lp.termEquals("", 1), "Empty member input ID");
				if(lp.termEquals("assigned", 3)){
					final long id=unsignedField(lp, 2);
					require(id<=Integer.MAX_VALUE, "Member family ID overflow");
					final int rank=ranks.get((int)id);
					require(rank>=0, "Assigned row references a family absent from codec");
					observed[rank]++;
					assigned++;
				}else{
					require(lp.termEquals("-1", 2), "Unassigned/filtered member must use family ID -1");
					if(lp.termEquals("unassigned", 3)){
						unassigned++;
					}else if(lp.termEquals("filtered_length", 3)){
						filtered++;
					}else{
						throw new IllegalArgumentException("Unknown member status");
					}
				}
			}
		}finally{require(!input.close(), "Member table read failed");}
		require(ordinal==expected, "Member table ended before declared population");
		for(int i=0; i<observed.length; i++){
			require(observed[i]==loaded.families[i].members, "Member count differs for family "+loaded.families[i].familyId);
		}
		return new long[]{assigned, unassigned, filtered};
	}

	/** Strict decimal parsing over native field bounds rejects overflow without per-row Strings. */
	private static long unsignedField(final LineParser1 lp, final int term){
		lp.setBounds(term);
		long value=0;
		require(lp.b()>lp.a(), "Empty integer field");
		for(int i=lp.a(); i<lp.b(); i++){
			final int digit=lp.line()[i]-'0';
			require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10, "Invalid or overflowing decimal field");
			value=value*10+digit;
		}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Metadata and binding     ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the producer key/value schema, rejecting duplicate or malformed records. */
	private static HashMap<String, String> stats(final String path){
		final HashMap<String, String> values=new HashMap<String, String>();
		final ByteFile input=ByteFile.makeByteFile1(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(input.nextLine(), STATS_HEADER), "Wrong clustering stats header");
			for(byte[] line; (line=input.nextLine())!=null; ){
				lp.set(line);
				require(lp.terms()==2, "Malformed stats row");
				final String key=lp.parseString(0), value=lp.parseString(1);
				require(!key.isEmpty() && !value.isEmpty() && values.put(key, value)==null, "Empty or duplicate stats key/value");
			}
		}finally{require(!input.close(), "Stats read failed");}
		return values;
	}

	/** Requires explicit unique audit options; there are no inferred scientific defaults. */
	private static HashMap<String, String> options(final String[] args){
		final HashMap<String, String> o=new HashMap<String, String>();
		for(String arg:args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Options require key=value");
			final String key=arg.substring(0, split);
			require(OPTIONS.contains("|"+key+"|"), "Unknown audit option: "+key);
			require(o.put(key, arg.substring(split+1))==null, "Duplicate audit option: "+key);
		}
		return o;
	}

	/** Require a regular input file and apply the selected identity policy.
	 * SHA80 mode measures stored bytes against the required pin. Structural mode
	 * returns an optional supplied pin unchanged, or unmeasured when absent; it
	 * does not measure content or check size/modification-time stability. */
	private static String verify(final HashMap<String, String> o, final String key, final VerificationMode verification){
		final String path=required(o, key), pin=o.get(key+"sha80");
		require(Files.isRegularFile(Paths.get(path)), "Missing regular "+key+" file");
		if(pin==null){
			require(verification==VerificationMode.STRUCTURAL, "Missing "+key+"sha80");
			return "unmeasured";
		}
		require(pin.matches("[0-9a-f]{20}"), "Invalid "+key+" sha80");
		if(verification==VerificationMode.SHA80){
			require(DigestSuffix.file(path).equals(pin), "Changed or incorrectly pinned "+key+" file");
		}
		return pin;
	}

	/** Remove and parse verification, case-insensitively; default to SHA80. */
	private static VerificationMode verification(final HashMap<String, String> o){
		final String value=o.containsKey("verification") ? o.remove("verification") : "sha80";
		if("sha80".equalsIgnoreCase(value)){return VerificationMode.SHA80;}
		if("structural".equalsIgnoreCase(value)){return VerificationMode.STRUCTURAL;}
		throw new IllegalArgumentException("verification must be sha80 or structural");
	}

	/** Missing metadata must fail instead of becoming a zero or an inferred default. */
	private static String required(final HashMap<String, String> o, final String key){
		final String value=o.get(key);
		require(value!=null && !value.isEmpty(), "Missing "+key);
		return value;
	}
	/** Decimal counts are nonnegative, with Long.parseLong providing overflow rejection. */
	private static long number(final HashMap<String, String> o, final String key){
		final String text=required(o, key);
		require(text.matches("[0-9]+"), "Expected nonnegative decimal for "+key);
		return Long.parseLong(text);
	}
	/** Read a nonnegative decimal count that must fit in a signed int. */
	private static int integer(final HashMap<String, String> o, final String key){
		return Math.toIntExact(number(o, key));
	}
	/** Reject malformed or inconsistent audit input independently of assertions. */
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Three files checked before and after semantic validation, in receipt order. */
	private static final String[] FILES={"codec", "members", "stats"};
	/** Exact four-column assignment schema consumed by members(). */
	private static final byte[] MEMBER_HEADER="input_index\tinput_id\tfamily_id\tstatus".getBytes(StandardCharsets.US_ASCII);
	/** Exact producer summary schema consumed by stats(). */
	private static final byte[] STATS_HEADER="key\tvalue".getBytes(StandardCharsets.US_ASCII);
	/** Delimited case-sensitive option names; delimiters prevent substring matches. */
	private static final String OPTIONS="|codec|codecsha80|members|memberssha80|stats|statssha80|k|aabits|width|sketchseed|positionmode|partitions|overlap|blacklistsha80|inputs|minlen|maxalignlen|memberranking|verification|";
}
