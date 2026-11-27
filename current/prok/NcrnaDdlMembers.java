package prok;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ListNum;

/** Loads source-bound cobalamin construction members for DDL profiles.
 * Joins the original clan routing ledger, builder assignments and actual FASTA
 * by complete headers. Model IDs bind to the supplied current library order;
 * header descriptions do not substitute for the published member census.
 * Original small/orphan records are checked but never acquire invented votes.
 *
 * One owner loads once, before ranking workers start. Only retained sequence
 * arrays survive loading; the full-ID maps are then discarded. Returned arrays
 * are borrowed read-only by fresh native sketches, not copied per geometry.
 * Current subset loads accept explicit singleton-addition ledgers without
 * replacing original assignments. Historical loaders retain their RF00174 schema.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlMembers{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the complete frozen construction join once, retaining only kept sequences.
	 * Model names and reference sequences must have the current library's order.
	 * Opens and closes its own native input streams; malformed or incomplete joins
	 * throw rather than returning a partial profile population.
	 * @return Member arrays borrowed read-only by the experiment */
	public static NcrnaDdlMembers load(final String fasta, final String ledger,
			final String assignments, final String census, final String summary,
			final String[] names, final byte[][] references){
		return load(fasta, ledger, assignments, census, summary, names, references, null);
	}

	/** Optional provenance observer receives only actual kept voters during the
	 * same checked FASTA pass. Names and sequence arrays are borrowed read-only.
	 * A failing observer invalidates the load; the reader still drains safely.
	 */
	public static NcrnaDdlMembers load(final String fasta, final String ledger,
			final String assignments, final String census, final String summary,
			final String[] names, final byte[][] references, final KeptMemberSink sink){
		final NcrnaDdlMembers result=new NcrnaDdlMembers(names, references);
		final Loader loader=new Loader(result, sink);
		loader.summary(summary);
		loader.census(census, references);
		loader.ledger(ledger);
		loader.assignments(assignments);
		loader.sequences(fasta);
		return result;
	}

	/** Provenance callback; never receives excluded small/orphan input records. */
	public interface KeptMemberSink{
		void accept(String fullId, int model, byte[] bases);
	}

	/** New all-gene construction uses the actual builder assignment table, not
	 * historical clan routing, taxonomy, a previous census or old library members.
	 * Full headers bind every input record, including explicit nonvoting outcomes. */
	public static NcrnaDdlMembers loadRoster(final String fasta, final String assignments,
			final String[] names, final byte[][] references, final int expectedInput, final int expectedKept){
		return loadRoster(fasta, assignments, names, references, expectedInput, expectedKept, false, null, null);
	}

	/** Loads the actual retained subset, then applies a fully observed holdout.
	 * Removed clusters never acquire reassigned votes. All original FASTA rows,
	 * including removed/small/orphan records, reach observe before any include. */
	public static NcrnaDdlMembers loadSubsetRoster(final String fasta, final String assignments,
			final String[] names, final byte[][] references, final int expectedInput, final RosterFilter filter){
		require(filter!=null, "Subset profiling requires an explicit source-bound holdout policy");
		return loadSubsetRoster(fasta, assignments, null, names, references, expectedInput, filter);
	}

	/** Explicit singleton edges supplement original assignments without reassignment.
	 * Input dispositions count unique originals; added/excluded votes count edges.
	 * The source-bound filter observes and tests each original at most once. */
	public static NcrnaDdlMembers loadSubsetRoster(final String fasta, final String assignments,
			final String additions, final String[] names, final byte[][] references,
			final int expectedInput, final RosterFilter filter){
		require(filter!=null, "Subset profiling requires an explicit source-bound holdout policy");
		return loadRoster(fasta, assignments, names, references, expectedInput, -1, true, filter, additions);
	}

	/** One explicit clan member per resource; no inference from biological lengths. */
	public static String memberFamily(final String[] names){
		require(names!=null && names.length>0 && names[0]!=null, "A member resource needs its actual ordered model names");
		for(final String family:new String[]{"RF00174", "RF01482", "RF01689"}){
			final String prefix="cobalamin_"+family+"_";
			if(!names[0].startsWith(prefix)){continue;}
			for(final String name:names){
				require(name!=null && name.startsWith(prefix), "A profile resource cannot mix clan member identities");
			}
			return family;
		}
		throw new IllegalArgumentException("Unsupported cobalamin member model namespace: "+names[0]);
	}

	/** Two phases prevent FASTA order from leaking a later panel gene into votes. */
	public interface RosterFilter{
		void observe(String fullId, byte[] bases);
		default void finishObservation(){}
		/** Only an explicit source-bound supplement may admit a historical ID. */
		default boolean boundHistoricalSource(String fullId){return false;}
		boolean include(String fullId, byte[] bases);
	}

	private static NcrnaDdlMembers loadRoster(final String fasta, final String assignments,
			final String[] names, final byte[][] references, final int expectedInput, final int expectedKept,
			final boolean subset, final RosterFilter filter, final String additions){
		return loadRoster(fasta, assignments, names, references, expectedInput, expectedKept, subset, filter, additions, null);
	}

	/** Retained incidence observer for a checked training-GC construction join.
	 * Original observations/exclusions are unchanged; singleton edges remain explicit. */
	public static NcrnaDdlMembers loadSubsetRoster(final String fasta, final String assignments,
			final String additions, final String[] names, final byte[][] references,
			final int expectedInput, final RosterFilter filter, final KeptMemberSink sink){
		require(filter!=null && sink!=null, "Training incidence capture needs explicit exclusions and a retained-edge observer");
		return loadRoster(fasta, assignments, names, references, expectedInput, -1, true, filter, additions, sink);
	}

	private static NcrnaDdlMembers loadRoster(final String fasta, final String assignments,
			final String[] names, final byte[][] references, final int expectedInput, final int expectedKept,
			final boolean subset, final RosterFilter filter, final String additions, final KeptMemberSink sink){
		require(expectedInput>0 && (subset || expectedKept>0 && expectedKept<=expectedInput) && !Shared.TRIM_READ_DESCRIPTION,
			"Declare the full original input and retained-voter counts and preserve complete FASTA descriptions");
		final NcrnaDdlMembers result=new NcrnaDdlMembers(names, references);
		final String family=memberFamily(names), prefix="cobalamin_"+family+"_";
		require(subset || family.equals(FAMILY), "Historical roster profiles retain their RF00174 binding");
		final ObjectIntMap<String> models=new ObjectIntMap<String>(String.class), pending=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> sourceStates=new ObjectIntMap<String>(String.class);
		for(int model=0; model<names.length; model++){
			final String id=RrnaResourceIO.first(names[model]);
			require(id.startsWith(prefix) && !models.contains(id), "Actual member models need unique family-qualified identities");
			models.put(id, model);
		}
		final int[] counts=new int[names.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(assignments)){
			RrnaResourceIO.header(input.reader, ASSIGNMENT_HEADER, assignments);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6 && p.termStartsWith(family+"_", 3), "Roster assignments must retain their actual member family and all six fields");
				final String id=p.parseString(0);
				final String source=RrnaResourceIO.first(id);
				require(!pending.contains(id) && !sourceStates.contains(source)
					&& (!subset || source.startsWith(family+"_input_") || filter!=null && filter.boundHistoricalSource(id)),
					"Each original construction source needs one actual builder assignment and full header");
				final int model;
				if(p.termEquals("kept", 5)){
					final String label=RrnaResourceIO.first(p.parseString(4));
					final int assigned=models.get(label);
					require(assigned>=0 || subset && label.startsWith(prefix),
						"A kept assignment must identify an actual current or explicitly removed model");
					model=assigned>=0 ? assigned : names.length+2;
					if(assigned>=0){
						counts[model]++;
						result.kept++;
					}else{result.removed++;}
				}else{
					require(p.termEquals("", 4), "Nonvoting assignments must not claim a current model");
					if(p.termEquals("small", 5)){
						model=names.length;
						result.small++;
					}else{
						require(p.termEquals("orphan", 5), "Unknown builder state cannot acquire model votes");
						model=names.length+1;
						result.orphan++;
					}
				}
				pending.put(id, model);
				sourceStates.put(source, model);
				result.inputMembers++;
			}
		}
		require(result.inputMembers==expectedInput && (subset || result.kept==expectedKept)
			&& (long)result.kept+result.small+result.orphan+result.removed==result.inputMembers,
			"The actual assignment populations must conserve all declared input/voting records");
		final ObjectIntMap<String> extraModels=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> extraLengths=new ObjectIntMap<String>(String.class);
		if(additions!=null){
			readAdditions(additions, family, models, sourceStates, counts, extraModels, extraLengths);
			result.added=extraModels.size();
		}
		final String[][] headers=subset ? new String[names.length][] : null;
		for(int model=0; model<names.length; model++){
			require(counts[model]>0,
				"Every current model must bind original assignments or an explicit singleton edge before holdout: "+names[model]);
			result.members[model]=new byte[counts[model]][];
			if(subset){headers[model]=new String[counts[model]];}
		}
		result.assigned=subset ? counts.clone() : null;
		Arrays.fill(counts, 0);
		final FileFormat format=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		require(format.fasta(), "Roster construction input must retain its native FASTA framing");
		int observed=0;
		Throwable deferred=null;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 0))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch.list){
					if(deferred!=null){continue;}
					try{
						final int model=pending.get(read.id);
						require(model>=0 && read.mate==null && read.bases!=null && read.length()>0 && pending.remove(read.id),
							"FASTA full headers must match every original assignment exactly once, without missing/extra/non-nucleotide records");
						for(final byte base:read.bases){require("ACGTRYSWKMBDHVNacgtryswkmbdhvn".indexOf(base)>=0, "Construction sequences retain DNA IUPAC symbols without lossy conversion");}
						if(filter!=null){filter.observe(read.id, read.bases);}
						if(model<names.length){
							require(counts[model]<result.members[model].length, "A model cannot acquire votes beyond its actual builder assignments");
							if(subset){headers[model][counts[model]]=read.id;}
							result.members[model][counts[model]++]=read.bases;
						}
						final int extra=extraModels.get(RrnaResourceIO.first(read.id));
						if(extra>=0){
							require(read.length()==extraLengths.get(RrnaResourceIO.first(read.id))
								&& Arrays.equals(read.bases, references[extra]) && counts[extra]==0,
								"A retained singleton must use its exact original sequence once, without replacing a cluster vote");
							headers[extra][0]=read.id;
							result.members[extra][counts[extra]++]=read.bases;
						}
						observed++;
					}catch(RuntimeException | Error invalid){deferred=invalid;}
				}
				input.stream.returnList(batch);
			}
			// Native close joins producers; drain first so a malformed singleton
			// cannot leave a producer blocked on its full output queue.
			if(deferred instanceof RuntimeException){throw (RuntimeException)deferred;}
			if(deferred instanceof Error){throw (Error)deferred;}
		}
		require(observed==expectedInput && pending.isEmpty(), "Every assigned original construction record must occur in the actual FASTA");
		if(filter!=null){filter.finishObservation();}
		final ObjectIntMap<String> decisions=new ObjectIntMap<String>(String.class);
		for(int model=0; model<names.length; model++){
			require(counts[model]==result.members[model].length, "FASTA voters must conserve each actual model's builder membership");
			if(filter!=null){
				final byte[][] group=result.members[model];
				int kept=0;
				for(int member=0; member<group.length; member++){
					final String header=headers[model][member];
					int include=decisions.get(header);
					if(include<0){
						include=filter.include(header, group[member]) ? 1 : 0;
						decisions.put(header, include);
					}
					if(include==1){
						if(sink!=null){sink.accept(header, model, group[member]);}
						group[kept++]=group[member];
					}
					else{
						if(sourceStates.get(RrnaResourceIO.first(header))==model){
							result.excluded++;
							result.kept--;
						}else{result.excludedAdded++;}
					}
				}
				result.members[model]=kept==group.length ? group : Arrays.copyOf(group, kept);
			}
		}
		assert((long)result.kept+result.small+result.orphan+result.removed+result.excluded==result.inputMembers)
			: "Retained, removed, excluded and nonvoting assignments must conserve the original input after holdout";
		int votes=0;
		for(final byte[][] group:result.members){votes=Math.addExact(votes, group.length);}
		require(votes==result.votingMembers(), "Retained vote edges must conserve original voters plus explicit singleton edges after holdout");
		return result;
	}

	/** Checks the complete producer ledger; only currently retained singleton rows vote. */
	private static void readAdditions(final String path, final String family,
			final ObjectIntMap<String> models, final ObjectIntMap<String> sources, final int[] counts,
			final ObjectIntMap<String> extras, final ObjectIntMap<String> lengths){
		final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, ADDITION_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Singleton additions require the producer's complete six-column ledger");
				final String source=p.parseString(0), label=p.parseString(1);
				require(p.termEquals(family, 2) && p.termEquals(family, 3)
					&& p.termEquals("EXACT_SINGLE_MEMBER_NO_REASSIGNMENT", 5) && p.parseInt(4)>0
					&& sources.contains(source) && seen.add(source)
					&& label.equals("cobalamin_"+family+"_extra_"+source),
					"Singleton additions must bind distinct actual original sources and the producer's exact member/model recipe");
				final int model=models.get(label);
				if(model<0){continue;}
				require(counts[model]==0, "An added singleton cannot overwrite a model's original assignments");
				counts[model]=1;
				extras.put(source, model);
				lengths.put(source, p.parseInt(4));
			}
		}
	}

	private NcrnaDdlMembers(final String[] names, final byte[][] references){
		require(names!=null && references!=null && names.length>0 && names.length==references.length,
			"Construction members require the current ordered model names and sequences");
		members=new byte[names.length][][];
		for(int i=0; i<names.length; i++){
			require(names[i]!=null && references[i]!=null && references[i].length>0, "Every current model requires a name and sequence");
		}
		modelNames=names.clone();
	}

	/** Borrowed read-only member arrays in construction FASTA order for this model. */
	public byte[][] members(final int model){return members[model];}
	public int models(){return members.length;}
	public int inputMembers(){return inputMembers;}
	public int keptMembers(){return kept;}
	public int smallMembers(){return small;}
	public int orphanMembers(){return orphan;}
	public int removedMembers(){return removed;}
	public int excludedMembers(){return excluded;}
	/** All retained model-vote incidences, including explicit singleton edges. */
	public int votingMembers(){return kept+added-excludedAdded;}
	public int addedMembers(){return added;}
	public int excludedAddedMembers(){return excludedAdded;}
	public int assignedMembers(final int model){return assigned==null ? members[model].length : assigned[model];}

	/*--------------------------------------------------------------*/
	/*----------------       Construction Join      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Loader{
		Loader(final NcrnaDdlMembers result_, final KeptMemberSink sink_){
			result=result_;
			sink=sink_;
			final int n=result.models();
			counts=new int[n];
			minimum=new int[n];
			maximum=new int[n];
			for(int i=0; i<n; i++){
				final String id=RrnaResourceIO.first(result.modelNames[i]);
				require(id.startsWith("cobalamin_RF00174_") && !models.contains(id), "Current RF00174 model IDs must be unique: "+id);
				models.put(id, i);
			}
		}

		void summary(final String path){
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, SUMMARY_HEADER, path);
				final byte[] row=in.reader.nextLine();
				require(row!=null && in.reader.nextLine()==null, "Bind the single RF00174 construction-summary row");
				p.set(row);
				require(p.terms()==17 && p.termEquals(FAMILY, 0), "Unexpected construction-summary family or schema");
				result.inputMembers=p.parseInt(1);
				result.kept=p.parseInt(2);
				result.small=p.parseInt(3);
				result.orphan=p.parseInt(4);
				require(result.kept>0 && result.small>=0 && result.orphan>=0
					&& (long)result.kept+result.small+result.orphan==result.inputMembers
					&& p.parseInt(5)==0 && p.parseInt(6)==result.kept && p.parseInt(7)==result.models()
					&& p.parseInt(10)==0, "Summary must conserve original members and nonempty models without singleton additions");
			}
		}

		void census(final String path, final byte[][] references){
			long total=0;
			int rows=0;
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, NcrnaMemberLengthBounds.HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==19, "Census requires its complete nineteen-field schema");
					final int model=models.get(p.parseString(0));
					require(model>=0 && p.termEquals(FAMILY, 1), "Census must identify exactly the current RF00174 model set");
					final int n=p.parseInt(3), lo=p.parseInt(4), hi=p.parseInt(8);
					require(result.members[model]==null && p.parseInt(2)==references[model].length && n>0 && lo>=17 && hi>=lo,
						"Census model lengths, counts and extrema must match the current library");
					result.members[model]=new byte[n][];
					minimum[model]=lo;
					maximum[model]=hi;
					total+=n;
					rows++;
				}
			}
			require(rows==result.models() && total==result.kept, "Per-model census must conserve all kept construction members");
		}

		/** Reconstructs exactly the two header forms emitted by NcrnaClanTrainingPrep.write. */
		void ledger(final String path){
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, LEDGER_HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==8, "Clan ledger must preserve every routing and original-header field");
					if(!p.termEquals(FAMILY, 3)){continue;}
					final String source=p.parseString(1), biological=p.parseString(2), raw=p.parseString(7);
					require((source.equals(FAMILY) || source.equals("RF01482") || source.equals("RF01689"))
						&& RrnaResourceIO.first(raw).equals(p.parseString(0)), "Routing rows must retain their original source locus ID");
					final String header=raw+(source.equals(FAMILY) ? " member_family="+biological : "")+" build_group="+FAMILY;
					for(int i=0; i<header.length(); i++){
						require(header.charAt(i)<128, "Native FASTA full-ID matching requires the original ASCII construction headers");
					}
					require(!pending.contains(header) && p.parseInt(6)>=17, "Routing ledger members must have unique complete headers and seed-length sequences");
					pending.put(header, UNASSIGNED);
					lengths.put(header, p.parseInt(6));
				}
			}
			require(pending.size()==result.inputMembers, "Routing ledger must cover every original member, including cross-labelled and orphan records");
		}

		void assignments(final String path){
			final int[] states=new int[3];
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, ASSIGNMENT_HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==6, "Builder assignment requires its full six-field schema");
					final String header=p.parseString(0), label=p.parseString(4);
					require(pending.get(header)==UNASSIGNED && p.termStartsWith(FAMILY+"_", 3),
						"Each exact routed member needs one original builder assignment: "+header);
					final int state=p.termEquals("kept", 5) ? 0 : p.termEquals("small", 5) ? 1 : p.termEquals("orphan", 5) ? 2 : -1;
					require(state>=0, "Unknown builder disposition cannot acquire a model vote");
					final int model;
					if(state==0){
						model=models.get(RrnaResourceIO.first(label));
						require(model>=0, "Kept assignment names a model absent from the current library: "+label);
						counts[model]++;
					}else{
						require(label.isEmpty(), "Small/orphan records must not claim a model");
						model=result.models()+state-1;
					}
					pending.put(header, model);
					states[state]++;
				}
			}
			require(states[0]==result.kept && states[1]==result.small && states[2]==result.orphan,
				"Builder dispositions must match the complete published member summary");
			for(int i=0; i<counts.length; i++){
				require(counts[i]==result.members[i].length, "Assigned members differ from the published per-model census: "+result.modelNames[i]);
			}
			Arrays.fill(counts, 0);
		}

		/** Drains after a bad member before closing: native close is not cancellation. */
		void sequences(final String path){
			require(!Shared.TRIM_READ_DESCRIPTION, "Member joins require complete FASTA headers; do not trim descriptions");
			final FileFormat format=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			require(format.fasta(), "Construction members must use native FASTA framing");
			final int[] observedMin=new int[result.models()], observedMax=new int[result.models()];
			Arrays.fill(observedMin, Integer.MAX_VALUE);
			int records=0;
			Throwable deferred=null;
			try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1))){
				in.stream.start();
				for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
					for(final Read read:batch.list){
						if(deferred!=null){continue;}
						try{
							final int model=pending.get(read.id);
							require(model>=0 && read.mate==null && read.bases!=null && read.length()==lengths.get(read.id),
								"FASTA full header or length differs from the routed assignment: "+read.id);
							for(final byte base:read.bases){require(AminoAcid.baseToNumber[base&255]>=0, "Construction core contains an undefined base: "+read.id);}
							require(pending.remove(read.id) && lengths.remove(read.id), "Repeated FASTA member cannot vote twice: "+read.id);
							if(model<result.models()){
								require(counts[model]<result.members[model].length, "FASTA model membership exceeds its published census");
								result.members[model][counts[model]++]=read.bases;
								observedMin[model]=Math.min(observedMin[model], read.length());
								observedMax[model]=Math.max(observedMax[model], read.length());
								if(sink!=null){sink.accept(read.id, model, read.bases);}
							}
							records++;
						}catch(RuntimeException | Error invalid){deferred=invalid;}
					}
					in.stream.returnList(batch);
				}
				if(deferred instanceof RuntimeException){throw (RuntimeException)deferred;}
				if(deferred instanceof Error){throw (Error)deferred;}
			}
			require(records==result.inputMembers && pending.isEmpty() && lengths.isEmpty(), "FASTA must contain every routed/assigned member exactly once");
			for(int i=0; i<counts.length; i++){
				require(counts[i]==result.members[i].length && observedMin[i]==minimum[i] && observedMax[i]==maximum[i],
					"Actual member counts and lengths differ from census for "+result.modelNames[i]);
			}
			assert(records==result.kept+result.small+result.orphan) : "Voting and nonvoting dispositions must conserve the complete construction population";
		}

		final NcrnaDdlMembers result;
		final KeptMemberSink sink;
		final int[] counts, minimum, maximum;
		final ObjectIntMap<String> models=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> pending=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private int inputMembers, kept, small, orphan;
	private final String[] modelNames;
	private int[] assigned;
	private int removed, excluded, added, excludedAdded;
	private final byte[][][] members;
	private static final int UNASSIGNED=-2;
	private static final String FAMILY="RF00174";
	private static final String LEDGER_HEADER="locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header";
	private static final String ASSIGNMENT_HEADER="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus";
	private static final String ADDITION_HEADER="locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction";
	private static final String SUMMARY_HEADER="emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max";
}
