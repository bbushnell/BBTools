package prok;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Fixed or learned sense-oriented CM input flanks on an exact consensus axis.
 * Sequence identity is checked directly, without a runtime payload digest.
 * Parsed primitive arrays are private and shared read-only between workers.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmInputPadding{

	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes declared constants on every model, without inventing training support. */
	public static void writeFixed(final String path, final int pad5, final int pad3,
			final String[] names, final byte[][] sequences){
		axis(names, sequences);
		require(path!=null && pad5>=0 && pad3>=0,
			"Fixed input padding requires an output and nonnegative sense-oriented flanks");
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(path)){
			output.writer.println(FIXED_FORMAT);
			output.writer.print(row.append("#fixed\t").append(pad5).tab().append(pad3).nl());
			output.writer.println(FIXED_SCHEMA);
			for(int m=0; m<names.length; m++){
				row.clear().append(m).tab().append(NcrnaSupport.recordName(names[m])).tab();
				for(final byte base:sequences[m]){row.append(upper(base));}
				output.writer.print(row.nl());
			}
		}
	}

	/** Writes one selected rule; insufficient support always retains 20/20. */
	public static void write(final String path, final String rule, final int minCalls, final int minGroups,
			final String[] names, final byte[][] sequences, final int[] pad5, final int[] pad3,
			final int[] finiteCalls, final int[] sourceGroups){
		axis(names, sequences);
		require(path!=null && knownRule(rule) && minCalls>0 && minGroups>0
			&& pad5!=null && pad3!=null && finiteCalls!=null && sourceGroups!=null
			&& pad5.length==names.length && pad3.length==names.length
			&& finiteCalls.length==names.length && sourceGroups.length==names.length,
			"Input-padding export needs one rule, positive support minima and complete parallel model arrays");
		for(int m=0; m<names.length; m++){
			support(pad5[m], pad3[m], finiteCalls[m], sourceGroups[m], minCalls, minGroups);
		}
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(path)){
			output.writer.println(FORMAT);
			output.writer.print(row.append("#policy\t").append(rule).tab().append(minCalls).tab().append(minGroups).nl());
			output.writer.println(SCHEMA);
			for(int m=0; m<names.length; m++){
				row.clear().append(m).tab().append(NcrnaSupport.recordName(names[m])).tab();
				for(final byte base:sequences[m]){row.append(upper(base));}
				final boolean learned=finiteCalls[m]>=minCalls && sourceGroups[m]>=minGroups;
				output.writer.print(row.tab().append(pad5[m]).tab().append(pad3[m]).tab()
					.append(finiteCalls[m]).tab().append(sourceGroups[m]).tab().append(learned).nl());
			}
		}
	}

	/** Combines a bound percentile fit with separately selected ends and an explicit base. */
	public static void writeSelected(final String path, final int base5, final int base3,
			final Data source, final boolean[] selected5, final boolean[] selected3){
		require(path!=null && source!=null && !source.fixed && !source.selected && "P98_PLUS10".equals(source.rule)
			&& base5>=0 && base3>=0 && selected5!=null && selected3!=null
			&& selected5.length==source.size() && selected3.length==source.size(),
			"Selected padding combines one sequence-bound percentile fit with declared base flanks and complete end masks");
		for(int m=0; m<source.size(); m++){
			selection(selected5[m] ? source.pad5(m) : base5, selected3[m] ? source.pad3(m) : base3,
				source.finiteCalls(m), source.sourceGroups(m), source.minCalls, source.minGroups,
				base5, base3, selected5[m], selected3[m]);
		}
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(path)){
			output.writer.println(SELECTED_FORMAT);
			output.writer.print(row.append("#selected\t").append(source.rule).tab().append(source.minCalls)
				.tab().append(source.minGroups).tab().append(base5).tab().append(base3).nl());
			output.writer.println(SELECTED_SCHEMA);
			for(int m=0; m<source.size(); m++){
				row.clear().append(m).tab().append(NcrnaSupport.recordName(source.names[m])).tab();
				for(final byte base:source.sequences[m]){row.append(upper(base));}
				output.writer.print(row.tab().append(selected5[m] ? source.pad5(m) : base5)
					.tab().append(selected3[m] ? source.pad3(m) : base3).tab().append(source.finiteCalls(m))
					.tab().append(source.sourceGroups(m)).tab().append(selected5[m]).tab().append(selected3[m]).nl());
			}
		}
	}

	/** Every row binds its ordered record name and complete consensus bases. */
	public static Data read(final String path, final String[] names, final byte[][] sequences){
		axis(names, sequences);
		require(path!=null, "Selected input-padding resources need an explicit path");
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Input-padding resources need an explicit learned or fixed format header: "+path);
			p.set(header);
			require(p.terms()==1, "Input-padding format headers have one exact field: "+path);
			if(p.termEquals(FIXED_FORMAT, 0)){return readFixed(input, p, names, sequences);}
			if(p.termEquals(SELECTED_FORMAT, 0)){return readSelected(input, p, names, sequences);}
			require(p.termEquals(FORMAT, 0), "Unknown input-padding resource format: "+path);
			final byte[] policy=input.reader.nextLine();
			require(policy!=null, "Input-padding resources must declare their rule and support minima");
			p.set(policy);
			require(p.terms()==4 && p.termEquals("#policy", 0), "Input-padding policy records have four explicit columns");
			final String rule=p.parseString(1);
			final int minCalls=nonnegative(p, 2), minGroups=nonnegative(p, 3);
			require(knownRule(rule) && minCalls>0 && minGroups>0,
				"A selected input-padding resource must declare a known rule and positive support minima");
			RrnaResourceIO.header(input.reader, SCHEMA, path);
			final Data result=new Data(rule, minCalls, minGroups, names, sequences);
			int rows=0;
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(rows<names.length && p.terms()==8 && nonnegative(p, 0)==rows
					&& p.termEquals(NcrnaSupport.recordName(names[rows]), 1),
					"Input-padding rows must preserve every current consensus name in its exact library order");
				require(sameSequence(p, 2, sequences[rows]),
					"Input-padding sequence differs from the current consensus; refit instead of transplanting pads for "+names[rows]);
				final int pad5=nonnegative(p, 3), pad3=nonnegative(p, 4);
				final int calls=nonnegative(p, 5), groups=nonnegative(p, 6);
				support(pad5, pad3, calls, groups, minCalls, minGroups);
				final boolean learned=calls>=minCalls && groups>=minGroups;
				require(p.termEquals(learned ? "true" : "false", 7),
					"The learned marker must agree with both finite-call and distinct-source support");
				result.pad5[rows]=pad5;
				result.pad3[rows]=pad3;
				result.finiteCalls[rows]=calls;
				result.sourceGroups[rows]=groups;
				rows++;
			}
			require(rows==names.length, "A truncated padding resource cannot silently change missing models to zero padding");
			return result;
		}
	}

	private static Data readSelected(final Input input, final LineParser1 p,
			final String[] names, final byte[][] sequences){
		final byte[] policy=input.reader.nextLine();
		require(policy!=null, "Selected padding must declare the source percentile policy and both base flanks");
		p.set(policy);
		require(p.terms()==6 && p.termEquals("#selected", 0) && p.termEquals("P98_PLUS10", 1),
			"Selected padding declares its percentile rule, two fit minima and two sense-oriented base flanks");
		final String rule=p.parseString(1);
		final int minCalls=nonnegative(p, 2), minGroups=nonnegative(p, 3);
		final int base5=nonnegative(p, 4), base3=nonnegative(p, 5);
		require(minCalls>0 && minGroups>0, "A selected resource retains the original positive percentile-fit support minima");
		RrnaResourceIO.header(input.reader, SELECTED_SCHEMA, "selected input padding");
		final Data result=new Data(rule, minCalls, minGroups, names, sequences, true, base5, base3);
		int rows=0;
		for(byte[] line; (line=input.reader.nextLine())!=null; rows++){
			p.set(line);
			require(rows<names.length && p.terms()==9 && nonnegative(p, 0)==rows
				&& p.termEquals(NcrnaSupport.recordName(names[rows]), 1) && sameSequence(p, 2, sequences[rows]),
				"Selected padding binds every name and consensus once on the current ordered model axis");
			final int five=nonnegative(p, 3), three=nonnegative(p, 4);
			final int calls=nonnegative(p, 5), groups=nonnegative(p, 6);
			final boolean use5=flag(p, 7), use3=flag(p, 8);
			selection(five, three, calls, groups, minCalls, minGroups, base5, base3, use5, use3);
			result.pad5[rows]=five;
			result.pad3[rows]=three;
			result.finiteCalls[rows]=calls;
			result.sourceGroups[rows]=groups;
		}
		require(rows==names.length, "A selected resource cannot omit unsupported models or silently invent zero padding");
		return result;
	}

	private static boolean flag(final LineParser1 p, final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Selected end masks are explicit true or false");
		return p.termEquals("true", column);
	}
	private static void selection(final int pad5, final int pad3, final int calls, final int groups,
			final int minCalls, final int minGroups, final int base5, final int base3,
			final boolean selected5, final boolean selected3){
		require(pad5>=0 && pad3>=0 && calls>=0 && groups>=0 && groups<=calls && (calls==0)==(groups==0),
			"Selected padding preserves real nonnegative percentile-fit support and valid source-group counts");
		require(!selected5 && !selected3 || calls>=minCalls && groups>=minGroups,
			"An unsupported percentile fit cannot become an apparent supported override of the declared base");
		require(selected5 ? pad5!=base5 : pad5==base5, "Unselected five-prime ends retain their declared base; unchanged pads are not overrides");
		require(selected3 ? pad3!=base3 : pad3==base3, "Unselected three-prime ends retain their declared base; unchanged pads are not overrides");
	}

	/** Fixed records bind only the model axis; no learned-support columns exist. */
	private static Data readFixed(final Input input, final LineParser1 p,
			final String[] names, final byte[][] sequences){
		assert(input!=null && p!=null && names.length==sequences.length) :
			"The fixed header and validated current consensus axis precede fixed-policy decoding";
		final byte[] policy=input.reader.nextLine();
		require(policy!=null, "Fixed input padding must declare both sense-oriented flanks");
		p.set(policy);
		require(p.terms()==3 && p.termEquals("#fixed", 0), "Fixed padding declares exactly #fixed, pad5 and pad3");
		final int pad5=nonnegative(p, 1), pad3=nonnegative(p, 2);
		RrnaResourceIO.header(input.reader, FIXED_SCHEMA, "fixed input padding");
		final Data result=new Data(names, sequences, pad5, pad3);
		int rows=0;
		for(byte[] line; (line=input.reader.nextLine())!=null;){
			p.set(line);
			require(rows<names.length && p.terms()==3 && nonnegative(p, 0)==rows
				&& p.termEquals(NcrnaSupport.recordName(names[rows]), 1),
				"Fixed padding must bind every consensus name once in current library order");
			require(sameSequence(p, 2, sequences[rows]),
				"Fixed-padding consensus differs from the loaded library: "+names[rows]);
			rows++;
		}
		require(rows==names.length, "Truncated fixed padding cannot leave some current models on another policy");
		return result;
	}

	/** Native Parse.parseInt deliberately omits overflow detection. Resource
	 * dimensions must not wrap a malformed value into a plausible flank/count. */
	private static int nonnegative(final LineParser1 p, final int field){
		assert(p!=null && field>=0 && field<p.terms()) : "The enclosing resource schema supplies this numeric column";
		if(p.setBounds(field)==0){throw new IllegalArgumentException("Missing padding resource integer; column="+field);}
		int value=0;
		for(int i=p.a(); i<p.b(); i++){
			final int digit=p.line()[i]-'0';
			if(digit<0 || digit>9 || value>(Integer.MAX_VALUE-digit)/10){
				throw new IllegalArgumentException("Padding resource needs a nonnegative decimal integer within int range; column="+field);
			}
			value=value*10+digit;
		}
		return value;
	}

	private static boolean sameSequence(final LineParser1 p, final int field, final byte[] expected){
		assert(expected!=null && expected.length>0) : "The validated current consensus supplies the sequence identity being checked";
		if(p.setBounds(field)!=expected.length){return false;}
		final byte[] line=p.line();
		final int from=p.a();
		for(int i=0; i<expected.length; i++){
			if(upper(line[from+i])!=upper(expected[i])){return false;}
		}
		return true;
	}

	private static byte upper(final byte base){return base>='a' && base<='z' ? (byte)(base-32) : base;}

	private static void axis(final String[] names, final byte[][] sequences){
		require(names!=null && sequences!=null && names.length>0 && names.length==sequences.length,
			"Input-padding resources need the complete parallel model-name and consensus-sequence axis");
		for(int m=0; m<names.length; m++){
			require(names[m]!=null && !names[m].isEmpty() && sequences[m]!=null && sequences[m].length>0,
				"Every padding descriptor must bind an actual named consensus, including fallback models");
		}
	}

	private static void support(final int pad5, final int pad3, final int calls, final int groups,
			final int minCalls, final int minGroups){
		require(pad5>=0 && pad3>=0 && calls>=0 && groups>=0 && groups<=calls && (calls==0)==(groups==0),
			"Input flanks are nonnegative bases; finite support has positive source groups and empty support has none");
		require(calls>=minCalls && groups>=minGroups || pad5==20 && pad3==20,
			"Either support minimum failing must preserve the explicit 20-base fallback at both ends");
	}

	static boolean knownRule(final String rule){
		for(final String candidate:RULES){if(candidate.equals(rule)){return true;}}
		return false;
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Bound Data          ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Data{
		private Data(final String rule_, final int minCalls_, final int minGroups_, final String[] names_, final byte[][] sequences_){
			this(rule_, minCalls_, minGroups_, names_, sequences_, false, 20, 20);
		}
		private Data(final String rule_, final int minCalls_, final int minGroups_, final String[] names_,
				final byte[][] sequences_, final boolean selected_, final int base5_, final int base3_){
			final int models=names_.length;
			assert(models>0 && knownRule(rule_) && minCalls_>0 && minGroups_>0) :
				"Only a validated policy and nonempty current model axis may allocate serving padding arrays";
			rule=rule_;
			minCalls=minCalls_;
			minGroups=minGroups_;
			fixed=false;
			selected=selected_;
			base5=base5_;
			base3=base3_;
			names=names_;
			sequences=sequences_;
			modelCount=models;
			pad5=new int[models];
			pad3=new int[models];
			finiteCalls=new int[models];
			sourceGroups=new int[models];
		}
		private Data(final String[] names_, final byte[][] sequences_, final int five, final int three){
			final int models=names_.length;
			assert(models>0 && five>=0 && three>=0) :
				"Validated fixed flanks apply to every record on the complete current consensus axis";
			rule="FIXED";
			minCalls=minGroups=0;
			fixed=true;
			selected=false;
			names=names_;
			sequences=sequences_;
			modelCount=models;
			base5=five;
			base3=three;
			pad5=pad3=finiteCalls=sourceGroups=null;
		}
		public int size(){return modelCount;}
		public int pad5(final int model){
			assert(model>=0 && model<modelCount) : "Padding indexes the exact current model axis checked at resource load";
			return fixed ? base5 : pad5[model];
		}
		public int pad3(final int model){
			assert(model>=0 && model<modelCount) : "Padding indexes the exact current model axis checked at resource load";
			return fixed ? base3 : pad3[model];
		}
		public int finiteCalls(final int model){
			assert(model>=0 && model<modelCount) : "Support belongs to the current model axis; fixed policy has no fitted observations";
			return fixed ? 0 : finiteCalls[model];
		}
		public int sourceGroups(final int model){
			assert(model>=0 && model<modelCount) : "Source counts belong to the current model axis; fixed policy has no fitted observations";
			return fixed ? 0 : sourceGroups[model];
		}
		/** True only for declared constants, independent of learned support. */
		public final boolean fixed;
		/** A base plus selected end overrides, not a newly fitted percentile table. */
		public final boolean selected;
		public final int base5, base3;
		public final String rule;
		public final int minCalls, minGroups;
		private final int modelCount;
		private final int[] pad5, pad3, finiteCalls, sourceGroups;
		// Borrow the validated immutable library axis; selected export cannot transplant to another axis.
		private final String[] names;
		private final byte[][] sequences;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final String[] RULES={"P99_PLUS_POPULATION_SD", "P99_PLUS10", "P99_9", "P90_PLUS10", "P95_PLUS10", "P98_PLUS10"};
	private static final String FORMAT="#BBTools CM input padding v1";
	private static final String SCHEMA="model_ordinal0\tmodel\tconsensus_sequence\tpad5\tpad3\ttraining_finite\ttraining_groups\tlearned";
	private static final String FIXED_FORMAT="#BBTools CM input padding fixed v1";
	private static final String FIXED_SCHEMA="model_ordinal0\tmodel\tconsensus_sequence";
	private static final String SELECTED_FORMAT="#BBTools CM input padding selected v1";
	private static final String SELECTED_SCHEMA="model_ordinal0\tmodel\tconsensus_sequence\tpad5\tpad3\ttraining_finite\ttraining_groups\tselected5\tselected3";
}
