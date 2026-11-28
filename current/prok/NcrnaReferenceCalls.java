package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import parse.LineParser1;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaTextIO.Input;
import prok.NcrnaObservationData.Symbols;

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
		this(rules_, false);
	}

	/** Optional saved-GFF projection of native typed RNA; generic rules stay explicit. */
	public NcrnaReferenceCalls(final Rule[] rules_, final boolean typed_){
		require(rules_!=null && rules_.length>0, "Query grading requires an explicit family projection");
		rules=rules_.clone();
		typed=typed_;
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
		final Result result=result(input);
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
					if(typed){retain(result, parse(line, lineNumber, path));}
				}else if(parser.termEquals("rRNA", 2)){
					result.rrna++;
					if(typed){retain(result, parse(line, lineNumber, path));}
				}else{
					check(parser.termEquals("RNA", 2), "Unknown CallGenes feature type", path, lineNumber);
					retain(result, parse(line, lineNumber, path));
				}
			}
		}
		assert(result.features==result.cds+result.trna+result.rrna+result.included+result.excluded-result.projectedTyped)
			: "Every CallGenes feature must survive in a named query population counter";
		return result;
	}

	private static void retain(final Result result, final Call call){
		assert(result!=null && call!=null) : "Every projected native row retains its original line and population";
		result.calls.add(call);
		result.ruleCounts[call.ruleIndex]++;
		if(!call.nativeType.equals("RNA")){result.projectedTyped++;}
		if(call.included){result.included++;}
		else{result.excluded++;}
	}

	/** Creates the same explicit input population used by the saved-GFF reader. */
	public Result result(final String input){
		require(input!=null && !input.isEmpty(), "Live calls require the original manifest input identity");
		return new Result(input, rules.length+(typed ? TYPED_FAMILIES.length : 0));
	}
	/** Parses one borrowed native generic-RNA GFF row without changing population
	 * counters. Returns the same owned call as read(); non-RNA rows are rejected.
	 * Useful to consumers that handle CDS rows in the same streaming input pass.
	 */
	public Call parseReport(final byte[] line, final long lineNumber, final String path){
		parser.set(line);
		check(parser.terms()==9 && parser.termEquals("BBTools", 1) && parser.termEquals("RNA", 2),
			"Expected one native generic-RNA GFF row", path, lineNumber);
		return parse(line, lineNumber, path);
	}

	/** Projects a native final path without serializing an intermediate GFF. Uses
	 * Orf.appendGff's one-based bounds, first-token contig and display-family label.
	 * The result retains only generic RNA loci, just like read(); all types count.
	 */
	public void add(final Result result, final Orf orf){
		assert(result!=null && orf!=null) : "Live grading consumes actual final-path records";
		addReport(result, orf.type, orf.scafName, orf.start, orf.stop, orf.strand, orf.trnaModel, orf.ncrnaFamily);
	}

	/** Projects native report metadata without creating a temporary Orf. Bounds
	 * are zero-based inclusive and already forward-genomic, including reverse
	 * reports; this method never applies a strand flip. The input name follows
	 * Orf.appendGff's first-token/ASCII contract. A snapshot consumer may call this
	 * once per potential candidate to prepare immutable query loci, but those are
	 * NOT selected calls. In that use, Call.line is a candidate ordinal rather
	 * than an emitted GFF line; keep selected trial order separately.
	 */
	public void addReport(final Result result, final int type, final String name,
			final int start, final int stop, final int strand, final String model, final String family){
		assert(result!=null) : "Report projection needs an explicitly bound input population";
		require(!typed, "Typed projection currently supports saved native GFFs only; live projection remains generic-only");
		result.features++;
		if(type==ProkObject.CDS){
			result.cds++;
			return;
		}
		if(type==ProkObject.tRNA){
			result.trna++;
			return;
		}
		if(type==ProkObject.r16S || type==ProkObject.r23S || type==ProkObject.r5S || type==ProkObject.r18S){
			result.rrna++;
			return;
		}
		require(type==ProkObject.RNA && model!=null, "Live query must be a modeled generic RNA");
		requireAscii(model);
		final int ruleIndex=projectModel(model);
		final Rule rule=rules[ruleIndex];
		if(family!=null){
			require(rule.displayLabel.equals(NcrnaFamilyRegistry.displayName(family)),
				"Live family label conflicts with model projection: "+family);
		}
		final String contig=reportContig(name);
		require(start>=0
			&& stop>=start && stop<Integer.MAX_VALUE && (strand==0 || strand==1),
			"Live query must fit the saved-GFF coordinate and strand contract");
		final Call call=new Call(new Locus(contig, rule.family, start+1, stop+1,
			strand==0 ? '+' : '-'), model, result.features, ruleIndex, rule.include);
		result.calls.add(call);
		result.ruleCounts[ruleIndex]++;
		if(call.included){result.included++;}
		else{result.excluded++;}
	}

	/** Literal first-token identity emitted by Orf.appendGff, with its ASCII limit. */
	public static String reportContig(final String name){
		require(name!=null, "A live generic RNA query needs a literal contig name");
		int end=0;
		while(end<name.length() && name.charAt(end)!=' ' && name.charAt(end)!='\t'){end++;}
		final String contig=end==name.length() ? name : name.substring(0, end);
		requireAscii(contig);
		require(!contig.isEmpty() && !contig.equals("."), "Live query must fit the saved-GFF contig contract");
		return contig;
	}

	/** Whether this explicit query projection includes a scoring family. */
	public boolean scoresFamily(final String family){
		if(typed){for(final String name:TYPED_FAMILIES){if(name.equals(family)){return true;}}}
		for(final Rule rule:rules){if(rule.include && rule.family.equals(family)){return true;}}
		return false;
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
	/** Immutable explicit family/include rule for a named generic-RNA model. */
	public Rule modelRule(final String model){return rules[projectModel(model)];}
	/** Same whole-prefix rule as the String route, on a borrowed model token.
	 * No row String, model String or Call is allocated. This read-only method may
	 * share immutable rules; unlike parseReport(), it never touches the parser.
	 */
	public Rule modelRule(final byte[] model, final int from, final int to){
		require(model!=null && from>=0 && from<to && to<=model.length,
			"Byte model projection requires a nonempty bounded model token");
		Rule selected=null;
		for(final Rule rule:rules){
			if(!starts(model, from, to, rule.prefixBytes)){continue;}
			final int after=from+rule.prefixBytes.length;
			if(after<to && model[after]!='_' && model[after]!=' ' && model[after]!='\t'){continue;}
			require(selected==null, "Ambiguous whole-prefix family projection for a saved model token");
			selected=rule;
		}
		require(selected!=null, "Unknown generic RNA model; add an explicit include/exclude mapping for the saved token");
		return selected;
	}

	/** Parses borrowed row bytes, allocating only strings retained by the query row. */
	private Call parse(final byte[] line, final long lineNumber, final String path){
		final boolean generic=parser.termEquals("RNA", 2);
		final boolean trna=parser.termEquals("tRNA", 2);
		check(generic || typed && (trna || parser.termEquals("rRNA", 2)),
			"Typed query projection requires typed=t and a native RNA feature", path, lineNumber);
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
		final int identityStart=parser.a();
		int identityStop=identityStart;
		while(identityStop<limit && line[identityStop]!=',' && line[identityStop]!=';'){identityStop++;}
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
		check(modelStart<0 || modelStop>modelStart, generic ? "Missing model: query attribute" : "Empty model: query attribute", path, lineNumber);
		final String model=modelStart<0 ? "." : text(line, modelStart, modelStop);
		if(generic){
			check(modelStart>=0, "Missing model: query attribute", path, lineNumber);
			final int ruleIndex=projectModel(model);
			final Rule rule=rules[ruleIndex];
			if(familyStart>=0){
				check(familyStop-familyStart==rule.labelBytes.length && starts(line, familyStart, familyStop, rule.labelBytes),
					"Query family: label conflicts with model mapping", path, lineNumber);
			}
			final String reportedFamily=typed && familyStart>=0 ? text(line, familyStart, familyStop) : ".";
			return new Call(new Locus(contig, rule.family, (int)start, (int)stop, strand),
				model, lineNumber, ruleIndex, rule.include, "RNA", reportedFamily);
		}
		final String reportedFamily=familyStart<0 ? "." : text(line, familyStart, familyStop);
		final String identity=text(line, identityStart, identityStop);
		final int typeIndex=typedIndex(identity);
		check(typeIndex>=0 && trna==(typeIndex==0), "Unknown or conflicting native typed RNA identity", path, lineNumber);
		final String family=TYPED_FAMILIES[typeIndex];
		check(familyStart<0 || reportedFamily.equals(identity), "Typed family: label conflicts with native type token", path, lineNumber);
		// A model need not exist, but a known generic model cannot contradict the native type.
		for(final Rule rule:rules){
			if(modelStart>=0 && wholePrefix(model, rule.modelPrefix)){
				check(rule.family.equals(family), "Typed model conflicts with explicit generic family mapping", path, lineNumber);
			}
		}
		return new Call(new Locus(contig, family, (int)start, (int)stop, strand),
			model, lineNumber, rules.length+typeIndex, true, identity, reportedFamily);
	}

	/** Native type token, not model naming, determines the existing truth union. */
	private static int typedIndex(final String identity){
		if(identity.equals("tRNA")){return 0;}
		if(identity.equals("5S")){return 1;}
		if(identity.equals("16S") || identity.equals("18S")){return 2;}
		if(identity.equals("23S")){return 3;}
		return -1;
	}
	private static boolean wholePrefix(final String model, final String prefix){
		assert(model!=null && prefix!=null) : "Explicit model projection compares retained nonnull tokens";
		if(!model.startsWith(prefix)){return false;}
		final int end=prefix.length();
		return end==model.length() || model.charAt(end)=='_' || model.charAt(end)==' ' || model.charAt(end)=='\t';
	}
	/** A fresh axis for opt-in typed populations; callers cannot mutate the shared labels. */
	public static String[] typedFamilies(){return TYPED_FAMILIES.clone();}

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
			familyBytes=family.getBytes(StandardCharsets.UTF_8);
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
		/** Intern the canonical scoring family without allocating bytes per report. */
		int internFamily(final Symbols symbols){return symbols.intern(familyBytes, 0, familyBytes.length);}
		/** Corroborate a borrowed display label using the saved-reader byte contract. */
		boolean matchesLabel(final byte[] row, final int from, final int to){
			return to-from==labelBytes.length && starts(row, from, to, labelBytes);
		}
		private final byte[] prefixBytes, labelBytes, familyBytes;
	}

	/** Original generic RNA row and its explicit projection, retaining duplicates. */
	public static final class Call{
		private Call(final Locus locus_, final String model_, final long line_, final int ruleIndex_, final boolean included_){
			this(locus_, model_, line_, ruleIndex_, included_, "RNA", ".");
		}
		private Call(final Locus locus_, final String model_, final long line_, final int ruleIndex_, final boolean included_,
			final String nativeType_, final String reportedFamily_){
			locus=locus_;
			model=model_;
			line=line_;
			ruleIndex=ruleIndex_;
			included=included_;
			nativeType=nativeType_;
			reportedFamily=reportedFamily_;
		}
		public final Locus locus;
		public final String model;
		public final long line;
		public final int ruleIndex;
		public final boolean included;
		public final String nativeType, reportedFamily;
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
		public long features, cds, trna, rrna, included, excluded, projectedTyped;
	}

	private final Rule[] rules;
	private final boolean typed;
	private static final String[] TYPED_FAMILIES={"tRNA", "5S", "SSU", "LSU"};
	private final LineParser1 parser=new LineParser1('\t');
	private static final byte[] MODEL="model:".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FAMILY="family:".getBytes(StandardCharsets.US_ASCII);
}
