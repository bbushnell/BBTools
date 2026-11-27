package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import aligner.CovarianceModel;
import aligner.CovarianceModelStateCorridors;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** One versioned, optionally compressed resource for a frozen consensus axis.
 * Every model has an M record, including models that retain uniform fallback.
 * Only observed states of supported models have S records. Four signed16-bit
 * offsets adjust the raw glocal start-low/high and end-low/high anchors; their
 * individual values need not be ordered because their anchors differ.
 * Names, member identities, state dimensions and record order are checked
 * against the caller's supplied axis. No file hashes enter the runtime API.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmCorridorResource{

	/** Writes complete model descriptors and only actual observed-state bounds. */
	public static void write(final String path, final String[] names, final String[] members,
			final CovarianceModel[] models, final CovarianceModelStateCorridors[] profiles,
			final int[] retainedCalls, final int[] resolvedGroups){
		write(path, names, members, models, profiles, retainedCalls, resolvedGroups, null, 20);
	}

	/** Explicit separate-context format. Complete consensus bases and anchor
	 * padding bind the fit; raw translated anchors receive residuals before clipping.
	 * The original writer continues emitting exactly the version1 format.
	 */
	public static void writeContext(final String path, final String[] names, final String[] members,
			final CovarianceModel[] models, final CovarianceModelStateCorridors[] profiles,
			final int[] retainedCalls, final int[] resolvedGroups, final byte[][] sequences, final int anchorPadding){
		sequenceAxis(names, sequences);
		require(anchorPadding>=0, "A separate-context fit must declare its nonnegative per-end anchor padding");
		write(path, names, members, models, profiles, retainedCalls, resolvedGroups, sequences, anchorPadding);
	}

	private static void write(final String path, final String[] names, final String[] members,
			final CovarianceModel[] models, final CovarianceModelStateCorridors[] profiles,
			final int[] retainedCalls, final int[] resolvedGroups, final byte[][] sequences, final int anchorPadding){
		axis(names, members, models);
		require(path!=null && profiles!=null && profiles.length==names.length && retainedCalls!=null
			&& retainedCalls.length==names.length && resolvedGroups!=null && resolvedGroups.length==names.length,
			"A corridor resource needs one profile and retained support descriptor per frozen consensus");
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(path)){
			output.writer.println(sequences==null ? FORMAT : CONTEXT_FORMAT);
			if(sequences!=null){
				output.writer.print(row.append("#anchor_padding\t").append(anchorPadding).nl());
				output.writer.println(CONTEXT_PROJECTION);
			}
			output.writer.println(sequences==null ? SCHEMA : CONTEXT_SCHEMA);
			for(int m=0; m<names.length; m++){
				final CovarianceModelStateCorridors profile=profiles[m];
				final int states=models[m]==null ? 0 : models[m].states();
				support(retainedCalls[m], resolvedGroups[m], profile!=null);
				require(profile==null || profile.belongsTo(models[m]),
					"Written bounds must belong to their exact declared CM, not another model with equal dimensions");
				row.clear().append('M').tab().append(m).tab().append(names[m]).tab().append(members[m])
					.tab().append(states).tab().append(retainedCalls[m]).tab().append(resolvedGroups[m]);
				if(sequences!=null){
					row.tab();
					for(final byte base:sequences[m]){row.append(upper(base));}
				}
				output.writer.print(row.nl());
				if(profile==null){continue;}
				int observed=0;
				for(int state=0; state<states; state++){
					if(!profile.observed(state)){continue;}
					observed++;
					row.clear().append('S').tab().append(state);
					for(int bound=0; bound<4; bound++){row.tab().append(profile.offset(state, bound));}
					output.writer.print(row.nl());
				}
				require(observed>0, "A supported trajectory fit cannot export a profile without any visited state");
			}
		}
	}

	/** Loads the resource onto the caller's existing CM objects. Missing state
	 * records remain unobserved and therefore keep the declared uniform fallback.
	 * Unknown/repeated/reordered models or states fail instead of being ignored.
	 */
	public static Data read(final String path, final String[] names, final String[] members,
			final CovarianceModel[] models){
		return read(path, names, members, models, null);
	}

	/** Version2 requires the actual current sequences. Legacy reads retain their
	 * twenty-base same-window contract; callers must inspect unclippedAnchors.
	 */
	public static Data read(final String path, final String[] names, final String[] members,
			final CovarianceModel[] models, final byte[][] sequences){
		axis(names, members, models);
		require(path!=null, "Corridor loading needs the explicitly selected resource path");
		final LineParser1 p=new LineParser1('\t');
		final Data result;
		int current=-1, previousState=-1;
		short[] offsets=null;
		boolean[] observed=null;
		try(Input input=new Input(path)){
			final byte[] format=input.reader.nextLine();
			final boolean context=Arrays.equals(format, CONTEXT_FORMAT_BYTES);
			require(context || Arrays.equals(format, FORMAT_BYTES), "Unknown state-corridor format cannot supply a guessed projection contract");
			int anchorPadding=20;
			if(context){
				sequenceAxis(names, sequences);
				final byte[] policy=input.reader.nextLine();
				require(policy!=null, "A separate-context corridor resource must declare its anchor padding");
				p.set(policy);
				require(p.terms()==2 && p.termEquals("#anchor_padding", 0), "Anchor-context policy has an explicit name and one base count");
				anchorPadding=p.parseInt(1);
				require(anchorPadding>=0, "Anchor padding cannot be negative or inferred from a search-input resource");
				RrnaResourceIO.header(input.reader, CONTEXT_PROJECTION, path);
			}
			result=new Data(names.length, anchorPadding, context);
			RrnaResourceIO.header(input.reader, context ? CONTEXT_SCHEMA : SCHEMA, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("M", 0)){
					if(current>=0){store(result, current, models[current], offsets, observed);}
					current++;
					require(current<names.length && p.terms()==(context ? 8 : 7) && p.parseInt(1)==current
						&& p.termEquals(names[current], 2) && p.termEquals(members[current], 3),
						"Corridor descriptors must conserve the exact ordered consensus names and member identities");
					if(context){
						require(sameSequence(p, sequences[current]), "A context corridor cannot transplant fitted state offsets onto a changed consensus sequence");
					}
					final int states=models[current]==null ? 0 : models[current].states();
					require(p.parseInt(4)==states, "Corridor state dimensions must match the caller's declared CM");
					result.retainedCalls[current]=p.parseInt(5);
					result.resolvedGroups[current]=p.parseInt(6);
					support(result.retainedCalls[current], result.resolvedGroups[current], result.retainedCalls[current]>0);
					require(states>0 || result.retainedCalls[current]==0,
						"An unrefined member without a CM cannot acquire fitted state bounds");
					offsets=new short[Math.multiplyExact(states, 4)];
					observed=new boolean[states];
					previousState=-1;
				}else{
					require(p.termEquals("S", 0) && p.terms()==6 && current>=0 && result.retainedCalls[current]>0,
						"Only a supported, already declared consensus may contain state records");
					final int state=p.parseInt(1);
					require(state>previousState && state<observed.length,
						"Observed state records must be unique, increasing and inside their own CM");
					previousState=state;
					observed[state]=true;
					for(int bound=0; bound<4; bound++){
						final int value=p.parseInt(bound+2);
						require(value>=Short.MIN_VALUE && value<=Short.MAX_VALUE,
							"Corridor residual offsets must fit signed16-bit storage without truncation");
						offsets[state*4+bound]=(short)value;
					}
				}
			}
		}
		require(current==names.length-1, "A truncated corridor resource cannot silently omit consensus fallback descriptors");
		store(result, current, models[current], offsets, observed);
		return result;
	}

	private static void sequenceAxis(final String[] names, final byte[][] sequences){
		require(names!=null && sequences!=null && names.length==sequences.length,
			"A context corridor needs the complete consensus-sequence axis, including uniform fallback descriptors");
		for(final byte[] sequence:sequences){
			require(sequence!=null && sequence.length>0, "Every context descriptor binds an actual nonempty consensus");
			for(final byte base:sequence){require(base>32 && base<127, "Consensus bases must fit one nonempty ASCII TSV field");}
		}
	}

	private static boolean sameSequence(final LineParser1 p, final byte[] sequence){
		assert(sequence!=null && sequence.length>0) : "The current model axis supplies the exact consensus being compared";
		if(p.setBounds(7)!=sequence.length){return false;}
		final int from=p.a();
		for(int i=0; i<sequence.length; i++){
			if(upper(p.line()[from+i])!=upper(sequence[i])){return false;}
		}
		return true;
	}

	private static byte upper(final byte base){return base>='a' && base<='z' ? (byte)(base-32) : base;}

	private static void store(final Data data, final int ordinal, final CovarianceModel model,
			final short[] offsets, final boolean[] observed){
		assert(ordinal>=0 && ordinal<data.profiles.length && observed!=null && offsets!=null) :
			"Only a fully read model descriptor can be finalized onto the resource's fixed axis";
		if(data.retainedCalls[ordinal]==0){return;}
		boolean any=false;
		for(final boolean value:observed){any|=value;}
		require(any && model!=null, "Supported corridor descriptors require at least one observed state on their actual CM");
		data.profiles[ordinal]=new CovarianceModelStateCorridors(model, offsets, observed);
	}

	private static void axis(final String[] names, final String[] members, final CovarianceModel[] models){
		require(names!=null && names.length>0 && members!=null && models!=null
			&& names.length==members.length && names.length==models.length,
			"Corridor resources require complete parallel consensus, member and CM axes");
		for(int i=0; i<names.length; i++){
			require(names[i]!=null && !names[i].isEmpty() && members[i]!=null && !members[i].isEmpty()
				&& names[i].indexOf('\t')<0 && names[i].indexOf('\n')<0
				&& members[i].indexOf('\t')<0 && members[i].indexOf('\n')<0,
				"Consensus/member identifiers must be nonempty single TSV fields");
			require(models[i]==null || members[i].equals(models[i].accession),
				"The supplied CM must belong to the declared biological member");
		}
	}

	private static void support(final int calls, final int groups, final boolean fitted){
		require(fitted ? calls>0 && groups>0 && groups<=calls : calls==0 && groups==0,
			"Uniform fallback has zero fitted support; learned bounds require positive retained calls and resolved groups");
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	public static final class Data{
		private Data(final int models, final int anchorPadding_, final boolean unclippedAnchors_){
			assert(models>0) : "The validated consensus axis supplies at least one model descriptor";
			profiles=new CovarianceModelStateCorridors[models];
			retainedCalls=new int[models];
			resolvedGroups=new int[models];
			anchorPadding=anchorPadding_;
			unclippedAnchors=unclippedAnchors_;
		}
		public final CovarianceModelStateCorridors[] profiles;
		public final int[] retainedCalls, resolvedGroups;
		public final int anchorPadding;
		public final boolean unclippedAnchors;
	}

	private static final String FORMAT="#BBTools CM state corridors v1";
	private static final String SCHEMA="#M: ordinal,name,member,states,retained_calls,resolved_groups; S: state,start_low,start_high,end_low,end_high";
	private static final String CONTEXT_FORMAT="#BBTools CM state corridors v2";
	private static final String CONTEXT_PROJECTION="#projection\tTRANSLATE_RESIDUAL_THEN_CLIP";
	private static final String CONTEXT_SCHEMA="#M: ordinal,name,member,states,retained_calls,resolved_groups,consensus_sequence; S: state,start_low,start_high,end_low,end_high";
	private static final byte[] FORMAT_BYTES=FORMAT.getBytes(StandardCharsets.US_ASCII);
	private static final byte[] CONTEXT_FORMAT_BYTES=CONTEXT_FORMAT.getBytes(StandardCharsets.US_ASCII);
}
