package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import parse.LineParser1;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaTextIO.Input;

/** Reads CallGenes query GFFs through an explicit model-to-scoring-family map.
 * CallGenes uses literal contig names and comma-delimited model:/family: fields,
 * not percent-decoded reference GFF3 attributes (Orf.appendGff). A present family
 * label must corroborate the mapped model. Older rows without a family label are
 * supported; unknown or multiply mapped generic RNA models fail loudly.
 *
 * Each instance owns a reusable byte parser and belongs to one thread. Results
 * retain every generic RNA row, including explicitly excluded families and exact
 * duplicate calls. No biological truth, model confidence or FASTA provenance is
 * inferred. The enclosing grader must establish the manifest's input/file binding.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceCalls{

	/** Copies the explicit rule array; immutable rules may be shared between readers. */
	public NcrnaReferenceCalls(final Rule[] rules_){
		require(rules_!=null && rules_.length>0, "Query grading requires an explicit family projection");
		rules=rules_.clone();
		for(int i=0; i<rules.length; i++){
			require(rules[i]!=null, "Query family projection contains a null rule");
			for(int j=0; j<i; j++){
				require(!rules[i].modelPrefix.equals(rules[j].modelPrefix),
					"Duplicate model prefix in family projection: "+rules[i].modelPrefix);
			}
		}
	}

	/** Reads one declared input's calls. Input identity is a caller-supplied binding,
	 * not something recoverable from CallGenes' bare contig names. Unknown feature
	 * types/sources fail rather than allowing a different producer's GFF as queries.
	 */
	public Result read(final String path, final String input){
		require(input!=null && !input.isEmpty(), "A query file needs exact manifest input identity");
		final Result result=new Result(input, rules.length);
		try(Input in=new Input(path)){
			long lineNumber=0;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				check(parser.terms()==9 && parser.termEquals("BBTools", 1),
					"Expected a nine-field CallGenes GFF row", path, lineNumber);
				result.features++;
				if(parser.termEquals("CDS", 2)){
					result.cds++;
				}else if(parser.termEquals("tRNA", 2)){
					result.trna++;
				}else if(parser.termEquals("rRNA", 2)){
					result.rrna++;
				}else{
					check(parser.termEquals("RNA", 2), "Unknown CallGenes feature type", path, lineNumber);
					final Call call=parse(line, lineNumber, path);
					result.calls.add(call);
					result.ruleCounts[call.ruleIndex]++;
					if(call.included){result.included++;}
					else{result.excluded++;}
				}
			}
		}
		assert(result.features==result.cds+result.trna+result.rrna+result.included+result.excluded)
			: "Every CallGenes feature must survive in a named query population counter";
		return result;
	}

	/** Creates the same explicit input population used by the saved-GFF reader. */
	public Result result(final String input){
		require(input!=null && !input.isEmpty(), "Live calls require the original manifest input identity");
		return new Result(input, rules.length);
	}

	/** Projects a native final path without serializing an intermediate GFF. Uses
	 * Orf.appendGff's one-based bounds, first-token contig and display-family label.
	 * The result retains only generic RNA loci, just like read(); all types count.
	 */
	public void add(final Result result, final Orf orf){
		assert(result!=null && orf!=null) : "Live grading consumes actual final-path records";
		result.features++;
		if(orf.type==ProkObject.CDS){result.cds++;return;}
		if(orf.type==ProkObject.tRNA){result.trna++;return;}
		if(orf.type==ProkObject.r16S || orf.type==ProkObject.r23S
				|| orf.type==ProkObject.r5S || orf.type==ProkObject.r18S){result.rrna++;return;}
		require(orf.type==ProkObject.RNA && orf.trnaModel!=null, "Live query must be a modeled generic RNA");
		requireAscii(orf.trnaModel);
		final int ruleIndex=projectModel(orf.trnaModel);
		final Rule rule=rules[ruleIndex];
		if(orf.ncrnaFamily!=null){
			require(rule.displayLabel.equals(NcrnaFamilyRegistry.displayName(orf.ncrnaFamily)),
				"Live family label conflicts with model projection: "+orf.ncrnaFamily);
		}
		final String name=orf.scafName;
		require(name!=null, "A live generic RNA query needs a literal contig name");
		int end=0;
		while(end<name.length() && name.charAt(end)!=' ' && name.charAt(end)!='\t'){end++;}
		final String contig=end==name.length() ? name : name.substring(0, end);
		requireAscii(contig);
		require(!contig.isEmpty() && !contig.equals(".") && orf.start>=0
			&& orf.stop>=orf.start && orf.stop<Integer.MAX_VALUE && (orf.strand==0 || orf.strand==1),
			"Live query must fit the saved-GFF coordinate and strand contract");
		final Call call=new Call(new Locus(contig, rule.family, orf.start+1, orf.stop+1,
			orf.strand==0 ? '+' : '-'), orf.trnaModel, result.features, ruleIndex, rule.include);
		result.calls.add(call);
		result.ruleCounts[ruleIndex]++;
		if(call.included){result.included++;}
		else{result.excluded++;}
	}

	/** Orf.appendGff uses ByteBuilder's byte-oriented char/string append. Reject
	 * non-ASCII live identifiers until that writer has a UTF-8 contract; preserving
	 * Unicode only here would silently disagree with its saved-GFF counterpart.
	 */
	private static void requireAscii(final String text){
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)>127){
				throw new IllegalArgumentException("Live projection requires ASCII identifiers because Orf.appendGff is byte-oriented: "+text);
			}
		}
	}

	/** Returns the unique rule for a model, using the saved/live whole-prefix contract.
	 * Configuration checks may use this before calling; no counters are changed.
	 */
	public int projectModel(final String model){
		require(model!=null && !model.isEmpty(), "Model projection requires a named consensus");
		int selected=-1;
		for(int i=0; i<rules.length; i++){
			final String prefix=rules[i].modelPrefix;
			if(!model.startsWith(prefix)){continue;}
			final int after=prefix.length();
			if(after<model.length() && model.charAt(after)!='_' && model.charAt(after)!=' ' && model.charAt(after)!='\t'){continue;}
			require(selected<0, "Ambiguous model-prefix family projection: "+model);
			selected=i;
		}
		require(selected>=0, "Unknown generic RNA model; add an explicit include/exclude mapping: "+model);
		return selected;
	}

	/** Parses borrowed row bytes, allocating only strings retained by the query row. */
	private Call parse(final byte[] line, final long lineNumber, final String path){
		final long start=parser.parseLong(3), stop=parser.parseLong(4);
		check(start>=1 && stop>=start && stop<=Integer.MAX_VALUE,
			"Query coordinates must fit positive one-based closed matcher bounds", path, lineNumber);
		final char strand=parser.termEquals("+", 6) ? '+' : parser.termEquals("-", 6) ? '-' : '.';
		check(strand!='.', "Query matching requires explicit strand", path, lineNumber);
		parser.setBounds(0);
		final String contig=text(line, parser.a(), parser.b());
		check(!contig.isEmpty() && !contig.equals("."), "A generic RNA query needs a literal contig", path, lineNumber);
		parser.setBounds(8);
		final int limit=parser.b();
		int modelStart=-1, modelStop=-1, familyStart=-1, familyStop=-1;
		for(int from=parser.a(); from<limit; ){
			int end=from;
			while(end<limit && line[end]!=',' && line[end]!=';'){end++;}
			if(starts(line, from, end, MODEL)){
				check(modelStart<0, "Duplicate model: query attribute", path, lineNumber);
				modelStart=from+MODEL.length;
				modelStop=end;
			}else if(starts(line, from, end, FAMILY)){
				check(familyStart<0, "Duplicate family: query attribute", path, lineNumber);
				familyStart=from+FAMILY.length;
				familyStop=end;
			}
			from=end+1;
		}
		check(modelStart>=0 && modelStop>modelStart, "Missing model: query attribute", path, lineNumber);
		final String model=text(line, modelStart, modelStop);
		final int ruleIndex=projectModel(model);
		final Rule rule=rules[ruleIndex];
		if(familyStart>=0){
			check(familyStop-familyStart==rule.labelBytes.length && starts(line, familyStart, familyStop, rule.labelBytes),
				"Query family: label conflicts with model mapping", path, lineNumber);
		}
		return new Call(new Locus(contig, rule.family, (int)start, (int)stop, strand),
			model, lineNumber, ruleIndex, rule.include);
	}

	private static boolean starts(final byte[] line, final int from, final int end, final byte[] prefix){
		assert(from>=0 && from<=end && end<=line.length) : "Query tokens must stay within the validated attribute field";
		if(end-from<prefix.length){return false;}
		for(int i=0; i<prefix.length; i++){if(line[from+i]!=prefix[i]){return false;}}
		return true;
	}

	/** CallGenes values are literal, so a '%' never invokes GFF3 decoding here. */
	private static String text(final byte[] line, final int from, final int end){
		assert(from>=0 && from<=end && end<=line.length) : "Retained query text is bounded by its source field";
		final String text=new String(line, from, end-from, StandardCharsets.UTF_8);
		require(text.indexOf('\uFFFD')<0, "Invalid UTF-8/replacement character would change query identity");
		return text;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	/** Builds contextual diagnostics only on failure, never once per valid row. */
	private static void check(final boolean valid, final String reason, final String path, final long line){
		if(!valid){throw new IllegalArgumentException(reason+" at "+path+":"+line);}
	}

	/** Explicit projection: several distinct model prefixes may share one scoring
	 * family. Exclusion is explicit and counted; no unmatched-model fallback exists.
	 */
	public static final class Rule{
		public Rule(final String family_, final String modelPrefix_, final String displayLabel_, final boolean include_){
			require(family_!=null && !family_.isEmpty() && modelPrefix_!=null && !modelPrefix_.isEmpty()
				&& displayLabel_!=null && !displayLabel_.isEmpty(), "A query rule needs scoring family, model prefix and expected display label");
			family=family_;
			modelPrefix=modelPrefix_;
			displayLabel=displayLabel_;
			include=include_;
			prefixBytes=modelPrefix.getBytes(StandardCharsets.UTF_8);
			labelBytes=displayLabel.getBytes(StandardCharsets.UTF_8);
			require(text(prefixBytes, 0, prefixBytes.length).equals(modelPrefix)
				&& text(labelBytes, 0, labelBytes.length).equals(displayLabel), "Query rule text must survive UTF-8 encoding");
		}
		public final String family, modelPrefix, displayLabel;
		public final boolean include;
		private final byte[] prefixBytes, labelBytes;
	}

	/** Original generic RNA row and its explicit projection, retaining duplicates. */
	public static final class Call{
		private Call(final Locus locus_, final String model_, final long line_, final int ruleIndex_, final boolean included_){
			locus=locus_;
			model=model_;
			line=line_;
			ruleIndex=ruleIndex_;
			included=included_;
		}
		public final Locus locus;
		public final String model;
		public final long line;
		public final int ruleIndex;
		public final boolean included;
	}

	/** Caller-owned rows and population counts for one exact declared input. */
	public static final class Result{
		private Result(final String input_, final int rules){
			input=input_;
			ruleCounts=new long[rules];
		}
		public final String input;
		public final ArrayList<Call> calls=new ArrayList<Call>();
		public final long[] ruleCounts;
		public long features, cds, trna, rrna, included, excluded;
	}

	private final Rule[] rules;
	private final LineParser1 parser=new LineParser1('\t');
	private static final byte[] MODEL="model:".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FAMILY="family:".getBytes(StandardCharsets.US_ASCII);
}
