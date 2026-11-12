package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Encodes and reads explicitly adjudicated ncRNA reference rows, without choosing
 * biological dispositions. One instance owns reusable parsing/encoding buffers and
 * must stay on one thread. Returned rows own retained strings and borrow no buffers.
 * Input identity comes from the per-file manifest and must also match the row.
 * Complete optional CM attributes include cm_models plus score/GA/percent/basis/
 * display precision. Partial prototypes lacking cm_models are rejected, never
 * interpreted from a protected reference_model assignment. Rows without any CM
 * attributes remain supported unchanged. Only scored/dont_care rows enter this format;
 * context and removed evidence belong
 * in separate products. This codec does not approve a rule or create a reference set.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceGff{

	/** Reads one canonical per-input file, including an empty reference. Requires
	 * the version/input directives written by write(); rejects repeated locus IDs.
	 * The returned list and its rows are caller-owned. Close failures propagate.
	 * No contig-length or biological-disposition validation is inferred here.
	 */
	public ArrayList<Row> read(final String path, final String expectedInput){
		require(expectedInput!=null && !expectedInput.isEmpty(), "The reference manifest must identify its exact input");
		final ArrayList<Row> rows=new ArrayList<Row>();
		final ObjectMap<String, Row> ids=new ObjectMap<String, Row>(16, String.class, Row.class);
		try(Input in=new Input(path)){
			final byte[] version=in.reader.nextLine();
			require(version!=null && same(version, 0, version.length, VERSION), "Canonical reference GFF requires its version directive: "+path);
			final byte[] binding=in.reader.nextLine();
			require(binding!=null && binding.length>INPUT_PREFIX.length()
				&& same(binding, 0, INPUT_PREFIX.length(), INPUT_PREFIX), "Reference GFF requires an input directive, even when empty: "+path);
			require(expectedInput.equals(decoder.decode(binding, INPUT_PREFIX.length(), binding.length)),
				"Reference file input directive disagrees with the exact manifest: "+path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='#'){
					require(!same(line, 0, line.length, "##FASTA"), "Reference GFF contains rows only, not embedded FASTA: "+path);
					require(!(line.length>=INPUT_PREFIX.length() && same(line, 0, INPUT_PREFIX.length(), INPUT_PREFIX))
						&& !(line.length>=VERSION_PREFIX.length() && same(line, 0, VERSION_PREFIX.length(), VERSION_PREFIX)),
						"Repeated reference header could concatenate different input files: "+path);
					continue;
				}
				final Row row=parse(line, expectedInput);
				require(!ids.contains(row.id), "Duplicate reference locus ID would double-count a scored or masked locus: "+row.id);
				ids.put(row.id, row);
				rows.add(row);
			}
		}
		return rows;
	}

	/** Writes explicit dispositions only, never derives them. Validates the whole
	 * supplied list before opening a fresh output; existing files are not overwritten.
	 * Input binding remains present when rows is empty. Caller retains the list.
	 */
	public void write(final String path, final String input, final List<Row> rows){
		require(input!=null && !input.isEmpty() && rows!=null, "Writing a reference file requires an input identity and explicit rows");
		final ObjectMap<String, Row> ids=new ObjectMap<String, Row>(16, String.class, Row.class);
		for(final Row row:rows){
			require(row!=null && input.equals(row.input), "Every reference row must belong to the declared file input");
			require(!ids.contains(row.id), "Duplicate reference locus ID would double-count a scored or masked locus: "+row.id);
			ids.put(row.id, row);
		}
		output.clear().append(VERSION).append('\n').append(INPUT_PREFIX);
		escape(output, input).append('\n');
		try(Output out=new Output(path)){
			for(final Row row:rows){
				append(row, output);
				if(output.length()>=16384){
					out.writer.print(output);
					output.clear();
				}
			}
			if(output.length()>0){out.writer.print(output);}
		}
	}

	/** Rewrites one declared evidence pointer without decoding/re-encoding other
	 * attributes or scientific fields. Requires canonical file/input framing and
	 * explicit row input/status/evidence attributes. This is not biological or
	 * full schema validation. Line content is preserved; native text I/O writes LF
	 * terminators. Compression bytes are not promised to remain identical.
	 * The caller owns witness relocation and its provenance checks.
	 * @return Number of feature rows whose evidence pointer was replaced.
	 */
	public long rebase(final String source, final String destination, final String input,
		final String from, final String to){
		require(input!=null && !input.isEmpty() && from!=null && to!=null && !from.equals(to),
			"Rebasing requires an exact input identity and distinct declared evidence pointers");
		final byte[] inputBytes=escape(new ByteBuilder(), input).toBytes();
		final byte[] oldBytes=escape(new ByteBuilder(), from).toBytes();
		final byte[] newBytes=escape(new ByteBuilder(), to).toBytes();
		long changed=0;
		try(Input in=new Input(source)){
			final byte[] version=in.reader.nextLine(), binding=in.reader.nextLine();
			require(version!=null && same(version, 0, version.length, VERSION), "Rebasing requires a canonical reference GFF version");
			require(binding!=null && binding.length>INPUT_PREFIX.length()
				&& same(binding, 0, INPUT_PREFIX.length(), INPUT_PREFIX)
				&& input.equals(decoder.decode(binding, INPUT_PREFIX.length(), binding.length)),
				"Reference file input directive disagrees with the rebase manifest: "+source);
			try(Output out=new Output(destination)){
				out.writer.println(version);
				out.writer.println(binding);
				for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
					if(line.length==0 || line[0]=='#'){
						require(!same(line, 0, line.length, "##FASTA")
							&& !(line.length>=INPUT_PREFIX.length() && same(line, 0, INPUT_PREFIX.length(), INPUT_PREFIX))
							&& !(line.length>=VERSION_PREFIX.length() && same(line, 0, VERSION_PREFIX.length(), VERSION_PREFIX)),
							"Rebasing cannot concatenate reference headers or embedded FASTA");
						out.writer.println(line);
					}else{
						rebaseLine(line, inputBytes, oldBytes, newBytes);
						out.writer.print(output);
						changed++;
					}
				}
			}
		}
		return changed;
	}

	/** Finds the unique declared attribute in bytes; no transient field Strings. */
	private void rebaseLine(final byte[] line, final byte[] input, final byte[] from, final byte[] to){
		parser.set(line);
		require(parser.terms()==9 && parser.termEquals("ncRNA", 2), "Rebasing requires nine-field ncRNA reference rows");
		parser.setBounds(8);
		final int limit=parser.b();
		int first=-1, last=-1, inputCount=0, statusCount=0;
		for(int start=parser.a(); start<limit; ){
			int end=start;
			while(end<limit && line[end]!=';'){end++;}
			int equals=start;
			while(equals<end && line[equals]!='='){equals++;}
			require(equals>start && equals<end, "Reference attributes require name=value framing");
			final int value=equals+1;
			if(same(line, start, equals, "reference_evidence")){
				require(first<0 && sameBytes(line, value, end, from), "Every feature must have one exact declared old evidence pointer");
				first=value;
				last=end;
			}else if(same(line, start, equals, "input_id")){
				require(++inputCount==1 && sameBytes(line, value, end, input), "Feature input_id must match the exact rebase manifest");
			}else if(same(line, start, equals, "reference_status")){
				require(++statusCount==1 && (same(line, value, end, "scored") || same(line, value, end, "dont_care")),
					"Rebasing requires one explicit scored/dont_care status; raw source annotations are not references");
			}
			start=end+1;
		}
		require(first>=0 && inputCount==1 && statusCount==1, "Missing reference input/status/evidence attributes");
		assert(last>=first && last<=line.length) : "The evidence replacement must stay inside its validated attribute value";
		output.clear().append(line, 0, first).append(to).append(line, last, line.length-last).nl();
	}
	private static boolean sameBytes(final byte[] line, final int from, final int to, final byte[] expected){
		assert(from>=0 && from<=to && to<=line.length) : "Rebase comparisons stay inside parsed attribute bounds";
		if(to-from!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(line[from+i]!=expected[i]){return false;}}
		return true;
	}

	/** Creates a row from an explicit disposition. No score or length implies a mask.
	 * Model may be null when unresolved; evidence must identify its decision record.
	 * A mask requires a nonempty reason; a scored row must have no mask reason.
	 */
	public Row create(final String input, final Locus locus, final boolean mask,
		final String reason, final String model, final String evidence){
		return create(input, locus, mask, reason, model, evidence, null);
	}

	/** Adds optional CM measurements. Values follow the separate cm_models order;
	 * legacy rows without measurements retain their original encoding. */
	public Row create(final String input, final Locus locus, final boolean mask,
		final String reason, final String model, final String evidence, final Scores scores){
		require(input!=null && !input.isEmpty(), "Reference rows need exact manifest input identity");
		require(locus!=null && locus.seqid!=null && !locus.seqid.isEmpty()
			&& locus.family!=null && !locus.family.isEmpty(), "Reference rows need literal contig and resolved family");
		require(locus.strand=='+' || locus.strand=='-', "Reference matching requires an explicit strand");
		require(evidence!=null && !evidence.isEmpty(), "Reference disposition needs a decision/evidence pointer");
		require(mask ? reason!=null && !reason.isEmpty() : reason==null,
			"dont_care requires a reason; scored rows cannot carry a mask reason");
		require(model==null || !model.isEmpty(), "An unresolved model is null, not an empty identifier");
		utf8IfNeeded(evidence);
		if(reason!=null){utf8IfNeeded(reason);}
		if(model!=null){utf8IfNeeded(model);}
		identity.clear().append("ncrna|");
		component(identity, input).append('|');
		component(identity, locus.seqid).append('|');
		component(identity, locus.family).append('|');
		identity.append(locus.start).append('|').append(locus.stop).append('|').append(locus.strand);
		final String id=new String(identity.array, 0, identity.length(), StandardCharsets.UTF_8);
		return new Row(input, locus, id, mask, reason, model, evidence, scores);
	}

	/** Appends one GFF3 feature and newline to the caller-owned reusable output buffer.
	 * Values are escaped exactly once, including literal percent signs in identifiers.
	 */
	public void append(final Row row, final ByteBuilder out){
		require(row!=null && out!=null, "Encoding needs a reference row and an output buffer");
		escape(out, row.locus.seqid).append("\tBBTools\tncRNA\t").append(row.locus.start)
			.append('\t').append(row.locus.stop).append("\t.\t").append(row.locus.strand).append("\t.\tID=");
		escape(out, row.id).append(";input_id=");
		escape(out, row.input).append(";reference_family=");
		escape(out, row.locus.family).append(";reference_status=").append(row.mask ? "dont_care" : "scored");
		if(row.mask){escape(out.append(";dont_care_reason="), row.reason);}
		if(row.model!=null){escape(out.append(";reference_model="), row.model);}
		escape(out.append(";reference_evidence="), row.evidence);
		if(row.scores!=null){appendScores(row.scores, out);}
		out.append('\n');
	}

	/** Reads one feature line without its newline. Header/comment handling belongs to
	 * the file reader. Missing/duplicate status, bad escapes and mismatched input fail.
	 * Unknown attributes are ignored; required identity/status attributes are exact.
	 */
	public Row parse(final byte[] line, final String expectedInput){
		require(line!=null, "A reference feature line is required");
		parser.set(line);
		require(parser.terms()==9 && parser.termEquals("ncRNA", 2), "Expected a nine-field ncRNA reference GFF row");
		final String contig=field(0);
		final long start=parser.parseLong(3), stop=parser.parseLong(4);
		require(start>=1 && stop>=start && stop<=Integer.MAX_VALUE,
			"Reference coordinates must fit positive one-based closed matcher bounds: "+start+".."+stop);
		final char strand=parser.termEquals("+", 6) ? '+' : parser.termEquals("-", 6) ? '-' : '.';
		Arrays.fill(values, null);
		parser.setBounds(8);
		final int limit=parser.b();
		for(int from=parser.a(); from<limit; ){
			int end=from;
			while(end<limit && line[end]!=';'){end++;}
			int equal=from;
			while(equal<end && line[equal]!='='){equal++;}
			require(equal>from && equal<end, "Reference attributes require nonempty name=value fields");
			for(int i=0; i<KEYS.length; i++){
				if(same(line, from, equal, KEYS[i])){
					require(values[i]==null, "Duplicate reference attribute: "+KEYS[i]);
					values[i]=decoder.decode(line, equal+1, end);
					require(!values[i].isEmpty(), "Empty reference attribute: "+KEYS[i]);
					break;
				}
			}
			from=end+1;
		}
		require(expectedInput!=null && expectedInput.equals(values[1]), "Reference input_id disagrees with the exact file manifest");
		require("scored".equals(values[3]) || "dont_care".equals(values[3]),
			"Missing or unknown reference_status; never infer a mask from companion metadata");
		final Scores scores=readScores();
		final Row result=create(expectedInput, new Locus(contig, values[2], (int)start, (int)stop, strand),
			"dont_care".equals(values[3]), values[4], values[5], values[6], scores);
		require(result.id.equals(values[0]), "Reference ID disagrees with its input/contig/family/interval/strand tuple");
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------       CM Measurements        ----------------*/
	/*--------------------------------------------------------------*/

	/** Transfers ownership of both arrays to the returned measurement record.
	 * Bits are exact integer tenths, not internal unrounded scores. GA strings
	 * preserve the supplied decimal values. Long.MIN_VALUE means no current hit;
	 * in that case GA may be either known or '.'. The caller must not mutate the
	 * arrays afterward. No threshold or biological disposition is selected here. */
	public Scores scores(final String models, final long[] bitsTenths, final String[] gaBits){
		require(bitsTenths!=null && gaBits!=null && bitsTenths.length>0 && bitsTenths.length==gaBits.length,
			"CM score and GA arrays must have identical nonzero lengths");
		require(bitsTenths.length==modelTerms(models), "Every CM score must identify its corresponding measured model");
		require(bitsTenths[0]==Long.MIN_VALUE || models!=null && !models.equals("."), "Measured CM scores require explicit model identities");
		return new Scores(models==null ? "." : models, bitsTenths, gaBits);
	}

	/** Retains source numbers for exact recomputation; displayed percentages have six decimals. */
	private void appendScores(final Scores scores, final ByteBuilder out){
		escape(out.append(";cm_models="), scores.models);
		out.append(";cm_bits=");
		for(int i=0; i<scores.size(); i++){
			if(i>0){out.append(',');}
			if(scores.missing){out.append('.');}
			else{appendFixed(out, scores.bitsTenths[i], 10);}
		}
		out.append(";cm_ga=");
		for(int i=0; i<scores.size(); i++){
			if(i>0){out.append(',');}
			out.append(scores.gaBits[i]);
		}
		out.append(";cm_percent_ga=");
		for(int i=0; i<scores.size(); i++){
			if(i>0){out.append(',');}
			if(scores.missing){out.append('.');}
			else{appendFixed(out, scores.percentMillionths[i], 1000000);}
		}
		out.append(";cm_score_basis=").append(scores.missing ? "no_current_native_hit" : "printed_0.1_bit");
		out.append(";cm_percent_precision=6_decimal_places");
	}

	/** Optional score attributes are all-or-none; no missing number becomes zero. */
	private Scores readScores(){
		int present=0;
		for(int i=7; i<values.length; i++){if(values[i]!=null){present++;}}
		if(present==0){return null;}
		require(present==6 && values[11].equals("6_decimal_places"), "CM measurement attributes require all six fields and explicit percent precision");
		final boolean missing=values[10].equals("no_current_native_hit");
		require(missing || values[10].equals("printed_0.1_bit"), "Unknown CM score basis cannot imply unrounded precision");
		final int count=terms(values[7]);
		require(count==terms(values[8]) && count==terms(values[9]), "CM score/GA/percent lists must have equal lengths");
		final long[] bits=new long[count];
		final String[] ga=new String[count];
		int bitStart=0, gaStart=0, percentStart=0;
		for(int i=0; i<count; i++){
			final int bitEnd=termEnd(values[7], bitStart), gaEnd=termEnd(values[8], gaStart);
			final int percentEnd=termEnd(values[9], percentStart);
			ga[i]=values[8].substring(gaStart, gaEnd);
			if(missing){
				require(isMissing(values[7], bitStart, bitEnd) && isMissing(values[9], percentStart, percentEnd),
					"A hitless prior must keep score and percentage explicitly missing");
				bits[i]=Long.MIN_VALUE;
			}else{bits[i]=tenths(values[7], bitStart, bitEnd);}
			bitStart=bitEnd+1;
			gaStart=gaEnd+1;
			percentStart=percentEnd+1;
		}
		final Scores result=scores(values[12], bits, ga);
		if(!missing){
			percentStart=0;
			for(int i=0; i<count; i++){
				final int end=termEnd(values[9], percentStart);
				appendFixed(scoreCheck.clear(), result.percentMillionths[i], 1000000);
				require(equal(scoreCheck, values[9], percentStart, end), "CM percentage must match its stored score/GA at the declared display precision");
				percentStart=end+1;
			}
		}
		return result;
	}

	/** ByteBuilder's double formatter omits decimals for integers; this fixed-scale
	 * form preserves the score format and explicit percentage display precision. */
	static void appendFixed(final ByteBuilder out, final long value, final int scale){
		assert(value!=Long.MIN_VALUE && (scale==10 || scale==1000000)) : "Missing CM scores bypass numeric encoding; only tenths and percentage-millionths are serialized";
		if(value<0){out.append('-');}
		final long magnitude=value<0 ? -value : value;
		out.append(magnitude/scale).append('.');
		long remainder=magnitude%scale;
		for(int divisor=scale/10; divisor>0; divisor/=10){
			out.append((byte)('0'+remainder/divisor));
			remainder%=divisor;
		}
	}

	/** Computes score/GA*100 to six places by decimal long division, rounding halves
	 * away from zero. Remainders stay below the <=15-digit GA numerator, so scaling
	 * them by ten cannot overflow. No per-event decimal objects or binary rounding. */
	static long percentage(final long bits, final long ga, final int gaDecimals){
		assert(bits!=Long.MIN_VALUE && ga>0 && ga<1000000000000000L && gaDecimals>=0) : "Exact percent division consumes validated tenth-bit scores and positive <=15-digit GA decimals";
		final long magnitude=bits<0 ? -bits : bits;
		long quotient=magnitude/ga, remainder=magnitude%ga;
		for(int i=0; i<gaDecimals+7; i++){
			final long scaled=remainder*10, digit=scaled/ga;
			require(quotient<=(Long.MAX_VALUE-digit)/10, "CM percentage exceeds fixed-six-decimal storage");
			quotient=quotient*10+digit;
			remainder=scaled%ga;
		}
		if(remainder>=(ga+1)/2){
			require(quotient<Long.MAX_VALUE, "Rounded CM percentage exceeds fixed-six-decimal storage");
			quotient++;
		}
		return bits<0 ? -quotient : quotient;
	}
	private static int modelTerms(final String text){
		if(text==null){return 1;}
		boolean empty=true;
		int count=1, from=0;
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)==','){
				require(!empty, "CM measurements cannot identify an empty model-list term");
				require(i-from!=1 || text.charAt(from)!='.', "A missing model marker cannot share a measured-model list");
				empty=true;
				count++;
				from=i+1;
			}else{empty=false;}
		}
		require(!empty, "CM measurements cannot identify an empty model-list term");
		require(count==1 || text.length()-from!=1 || text.charAt(from)!='.', "A missing model marker cannot share a measured-model list");
		return count;
	}

	/** Parses one-decimal native score text exactly, without a floating-point boundary. */
	private static long tenths(final String text, final int from, final int to){
		require(to-from>=3 && text.charAt(to-2)=='.', "Printed CM bits require exactly one decimal place");
		long value=0;
		int i=from;
		final boolean negative=text.charAt(i)=='-';
		if(negative){i++;}
		require(i<to-2, "Printed CM score needs an integer part");
		for(; i<to; i++){
			if(i==to-2){continue;}
			final int digit=text.charAt(i)-'0';
			require(digit>=0 && digit<=9 && value<(1L<<49), "Invalid or unrepresentable printed CM score");
			value=value*10+digit;
		}
		return negative ? -value : value;
	}
	private static int terms(final String text){
		int count=1;
		if(text!=null){for(int i=0; i<text.length(); i++){if(text.charAt(i)==','){count++;}}}
		return count;
	}
	private static int termEnd(final String text, final int from){
		final int comma=text.indexOf(',', from);
		return comma<0 ? text.length() : comma;
	}
	private static boolean isMissing(final String text, final int from, final int to){return to-from==1 && text.charAt(from)=='.';}
	private static boolean equal(final ByteBuilder bytes, final String text, final int from, final int to){
		if(bytes.length()!=to-from){return false;}
		for(int i=0; i<bytes.length(); i++){if(bytes.array[i]!=text.charAt(from+i)){return false;}}
		return true;
	}

	/** Decodes a retained field using the existing strict GFF percent/UTF-8 reader. */
	private String field(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}

	/** Length prefixes count UTF-8 bytes, not UTF-16 characters. Common ASCII values
	 * append directly; Unicode takes the existing Java UTF-8 encoder only when needed.
	 */
	private static ByteBuilder component(final ByteBuilder out, final String text){
		assert(out!=null && text!=null) : "ID components come from validated input/contig/family identities";
		final byte[] unicode=utf8IfNeeded(text);
		out.append(unicode==null ? text.length() : unicode.length).append(':');
		return unicode==null ? out.append(text) : out.append(unicode);
	}

	/** Escapes to the unreserved ASCII subset, valid in both seqid and attributes. */
	private static ByteBuilder escape(final ByteBuilder out, final String text){
		assert(out!=null && text!=null) : "The encoder only receives required values or present optional attributes";
		final byte[] unicode=utf8IfNeeded(text);
		final int length=unicode==null ? text.length() : unicode.length;
		for(int i=0; i<length; i++){
			final int value=unicode==null ? text.charAt(i) : unicode[i]&255;
			if(value>='A' && value<='Z' || value>='a' && value<='z'
				|| value>='0' && value<='9' || value=='.' || value=='_' || value=='-'){
				out.append((byte)value);
			}else{
				out.append('%').append(HEX[value>>>4]).append(HEX[value&15]);
			}
		}
		return out;
	}

	/** Returns null for ASCII; rejects malformed surrogate text rather than replacing
	 * identity with '?'. ByteBuilder.append(String) is intentionally ASCII-only.
	 */
	private static byte[] utf8IfNeeded(final String text){
		assert(text!=null) : "Required reference identities and present optional values have already been validated";
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)>127){
				require(text.indexOf('\uFFFD')<0, "Unicode replacement characters are rejected by the strict GFF identity reader");
				final byte[] bytes=text.getBytes(StandardCharsets.UTF_8);
				require(new String(bytes, StandardCharsets.UTF_8).equals(text), "Invalid Unicode would change reference identity");
				return bytes;
			}
		}
		return null;
	}

	private static boolean same(final byte[] line, final int from, final int end, final String key){
		assert(from>=0 && from<=end && end<=line.length) : "Attribute key comparison stays inside validated ninth-field bounds";
		if(end-from!=key.length()){return false;}
		for(int i=from; i<end; i++){if(line[i]!=key.charAt(i-from)){return false;}}
		return true;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}

	/** Immutable explicit reference disposition; a locus is immutable too. */
	public static final class Row{
		private Row(final String input_, final Locus locus_, final String id_, final boolean mask_,
			final String reason_, final String model_, final String evidence_, final Scores scores_){
			input=input_;
			locus=locus_;
			id=id_;
			mask=mask_;
			reason=reason_;
			model=model_;
			evidence=evidence_;
			scores=scores_;
		}
		public final String input, id, reason, model, evidence;
		public final Locus locus;
		public final boolean mask;
		/** Optional measured evidence; null only for legacy rows without these attributes. */
		public final Scores scores;
	}

	/** Owned measurement arrays, aligned with cm_models order. This is independent of
	 * a protected prior's resolved reference_model. Missing evidence
	 * has an explicit sentinel and basis; it cannot acquire a numeric percentage. */
	public static final class Scores{
		private Scores(final String models_, final long[] bitsTenths_, final String[] gaBits_){
			models=models_;
			bitsTenths=bitsTenths_;
			gaBits=gaBits_;
			percentMillionths=new long[bitsTenths.length];
			missing=bitsTenths[0]==Long.MIN_VALUE;
			for(int i=0; i<bitsTenths.length; i++){
				require((bitsTenths[i]==Long.MIN_VALUE)==missing, "Measured hits and hitless-prior markers cannot share one measurement list");
				require(missing || bitsTenths[i]>-(1L<<50) && bitsTenths[i]<(1L<<50), "Printed scores must fit lossless tenth-bit storage");
				final String value=gaBits[i];
				require(value!=null && !value.isEmpty(), "Every measurement needs a GA value or explicit missing marker");
				if(value.equals(".")){
					require(missing, "A measured hit cannot have an unknown GA");
					percentMillionths[i]=Long.MIN_VALUE;
					continue;
				}
				int points=0, digits=0, decimals=0;
				long numerator=0;
				for(int j=0; j<value.length(); j++){
					final char c=value.charAt(j);
					if(c=='.'){points++;}
					else{
						require(c>='0' && c<='9', "GA must be a plain positive decimal, not encoded attributes or an exponent");
						digits++;
						require(digits<=15, "GA decimal exceeds the supported exact input precision");
						numerator=numerator*10+c-'0';
						if(points>0){decimals++;}
					}
				}
				require(points<=1 && digits>0, "GA must be one plain decimal value");
				require(numerator>0, "A captured GA must be positive");
				percentMillionths[i]=missing ? Long.MIN_VALUE : percentage(bitsTenths[i], numerator, decimals);
			}
		}
		public int size(){return bitsTenths.length;}
		public long bitsTenths(final int index){return bitsTenths[index];}
		public String gaBits(final int index){return gaBits[index];}
		public final boolean missing;
		public final String models;
		private final long[] bitsTenths;
		private final String[] gaBits;
		private final long[] percentMillionths;
	}

	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder identity=new ByteBuilder();
	private final ByteBuilder output=new ByteBuilder();
	private final ByteBuilder scoreCheck=new ByteBuilder();
	private final String[] values=new String[13];
	private static final String[] KEYS={"ID", "input_id", "reference_family", "reference_status",
		"dont_care_reason", "reference_model", "reference_evidence", "cm_bits", "cm_ga", "cm_percent_ga",
		"cm_score_basis", "cm_percent_precision", "cm_models"};
	private static final char[] HEX="0123456789ABCDEF".toCharArray();
	private static final String VERSION="##gff-version 3";
	private static final String VERSION_PREFIX="##gff-version ";
	private static final String INPUT_PREFIX="##reference-input-id ";
}
