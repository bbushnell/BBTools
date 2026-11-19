package prok;

import java.util.ArrayList;
import java.util.IdentityHashMap;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Source-qualified boundary evidence from actual accepted candidates to final output.
 * One Writer borrows the shared native output; each Worker owns its candidate map,
 * score copies and scratch for one calling thread. No genome sequence is retained.
 * Original coordinates stay sense-oriented while final coordinates are genomic.
 * Call finish only after final path selection/reporting, and clear in finally.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryTrace{

	/** Borrows an already-started writer. The CallGenes invocation closes it after joining workers. */
	public static final class Writer{
		public Writer(final ByteStreamWriter output_, final String target_){
			this(output_, target_, false);
		}
		/** Invocation-only output observes the refiner call and final selection, without score replay. */
		public Writer(final ByteStreamWriter output_, final String target_, final boolean invocationsOnly_){
			require(output_!=null && target_!=null && !target_.isEmpty(), "Trace output and an explicit target family are required");
			output=output_;
			target=target_;
			invocationsOnly=invocationsOnly_;
		}
		public boolean accepts(final NcrnaFamily family){
			final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family.name);
			return family.name.equals(target) || entry.group!=null && entry.group.id.equals(target);
		}
		/** Fail before workers if the explicit target has no loaded caller resources. */
		public void validate(final ArrayList<NcrnaFamily> families){
			validateTarget(target, families, invocationsOnly);
		}
		/** Input records include empty assemblies; ordinals match ordinary source-qualified GFF. */
		public void input(final int ordinal, final String path){
			require(ordinal>=0 && path!=null && path.indexOf('\t')<0 && path.indexOf('\n')<0, "Trace input needs its exact ordered assembly path");
			final ByteBuilder row=new ByteBuilder("I\t").append(ordinal).tab().append(path).nl();
			synchronized(output){output.print(row);}
			if(selectionOutput!=null){synchronized(selectionOutput){selectionOutput.print(row);}}
		}
		public Worker worker(final int ordinal, final int thread){return new Worker(this, ordinal, thread);}
		/** Optional independent selected-geometry stream; configure before workers start. */
		public void setSelectionOutput(final ByteStreamWriter selected){
			require(invocationsOnly && selected!=null && selectionOutput==null,
				"Selected geometry requires one explicit invocation-only trace output");
			selectionOutput=selected;
		}
		private final ByteStreamWriter output;
		private final String target;
		private final boolean invocationsOnly;
		private ByteStreamWriter selectionOutput;
	}

	/** Worker-confined owner. begin/finish/clear bracket one ordinary contig call. */
	public static final class Worker{
		private Worker(final Writer writer_, final int ordinal_, final int thread_){
			require(ordinal_>=0 && thread_>=0, "Trace IDs require actual input and worker ordinals");
			writer=writer_;
			ordinal=ordinal_;
			thread=thread_;
		}
		public boolean accepts(final NcrnaFamily family){return writer.accepts(family);}
		public Family family(final NcrnaFamily family){
			require(accepts(family), "Cannot attach an unselected family to this trace");
			validateFamily(family, writer.invocationsOnly);
			final Family adapter=new Family(this, family);
			families.add(adapter);
			return adapter;
		}
		public void begin(final String name, final int length_, final long readID_){
			require(fullName==null && pending.isEmpty() && records.isEmpty() && length_>=0 && readID_>=0,
				"Clear previous trace state before beginning a contig with its actual input read ID");
			fullName=name;
			contig=NcrnaSupport.recordName(name);
			require(!contig.isEmpty(), "Trace contig must have the GFF's nonempty first-token identifier");
			length=length_;
			readID=readID_;
		}
		/** Preserves selection geometry before the final in-place reporting transformation. */
		public void beforeReport(final ArrayList<Orf> selected, final int sameCap, final int oppositeCap){
			if(writer.selectionOutput!=null){
				require(fullName!=null, "Selected geometry must belong to the active input contig");
				NcrnaSelectionTrace.write(writer.selectionOutput, buffer, ordinal, readID, contig,
					writer.target, selected, sameCap, oppositeCap);
			}
		}
		/** Records the selected objects after their in-place strand/reporting transformations. */
		public void finish(final ArrayList<Orf> selected){
			require(fullName!=null && selected!=null, "Final trace links require an active contig and a completed path");
			for(int finalIndex=0; finalIndex<selected.size(); finalIndex++){
				final Orf orf=selected.get(finalIndex);
				final Record record=pending.get(orf);
				if(record==null){
					if(orf.ncrnaFamily!=null && matches(orf.ncrnaFamily, writer.target)){
						throw new IllegalStateException("Selected target call lacks its candidate-bound trace: "+contig+" "+orf.ncrnaFamily);
					}
					continue;
				}
				require(record.returned && record.complete && !record.selected, "Final selected calls must be distinct completed scavenger returns");
				record.selected=true;
				record.finalIndex=finalIndex;
				record.finalStart=orf.start;
				record.finalStop=orf.stop;
				require(record.family.family.name.equals(orf.ncrnaFamily) && record.model.equals(NcrnaSupport.recordName(orf.trnaModel))
					&& record.strand==orf.strand, "Final path cannot silently change a traced call's family/model/strand");
				if(record.known){
					final int start=record.strand==0 ? record.refinedStart : length-1-record.refinedStop;
					final int stop=record.strand==0 ? record.refinedStop : length-1-record.refinedStart;
					require(start==orf.start && stop==orf.stop, "Final genomic trace bounds must be the transformed observed refinement");
				}
			}
			int returned=0, retained=0, unknown=0, unavailable=0;
			for(final Record record:records){
				require(record.complete, "Each recorded boundary seam must finish before final-path evidence is written");
				if(record.known && !writer.invocationsOnly){
					final long replay=record.endpoints(0);
					require((int)(replay>>>32)==record.refinedStart && (int)replay==record.refinedStop,
						"Zero-margin trace reconstruction differs from observed inference for "+record.model+" on "+contig);
				}
				if(record.returned){returned++;}
				if(record.selected){retained++;}
				if(!record.known){unknown++;}
				if(!writer.invocationsOnly && record.eligible && record.count==0){unavailable++;}
				write(record);
			}
			buffer.clear().append("T\t").append(ordinal).tab().append(thread).tab().append(readID).tab().append(contig).tab().append(length).tab()
				.append(records.size()).tab().append(returned).tab().append(retained).tab().append(unknown).tab().append(unavailable).nl();
			synchronized(writer.output){writer.output.print(buffer);}
		}
		/** Always called in finally, including when a reader, alignment, path or writer fails. */
		public void clear(){
			pending.clear();
			records.clear();
			for(final Family family:families){
				family.current=null;
				family.windowOffset=0;
			}
			fullName=contig=null;
			length=0;
		}
		public boolean empty(){return fullName==null && pending.isEmpty() && records.isEmpty();}
		private Record add(final Family family, final Orf orf, final int model, final boolean known,
				final boolean trim, final boolean eligible, final int fixed, final float gc, final float acceptedIdentity, final float identity,
				final boolean quantum, final int min, final int max, final int windowStart, final int windowStop, final int alignedLength){
			require(fullName!=null && fullName.equals(orf.scafName) && length==orf.scaflen && orf.flipped()==0,
				"Candidate capture must precede strand flipping on its actual input contig");
			require(!pending.containsKey(orf) && model>=0 && model<family.family.modelNames.length,
				"A selected-model candidate may enter the trace only once");
			final Record record=new Record(family, orf, model, nextID, known, trim, eligible, fixed, gc, acceptedIdentity,
				identity, quantum, min, max, windowStart, windowStop, alignedLength);
			nextID=Math.addExact(nextID, 1);
			pending.put(orf, record);
			records.add(record);
			return record;
		}
		private void write(final Record r){
			buffer.clear().append("C\t").append(ordinal).tab().append(thread).tab().append(r.id).tab().append(readID).tab().append(contig).tab()
				.append(r.strand).tab().append(r.family.family.name).tab().append(r.model).tab().append(length).tab()
				.append(r.start).tab().append(r.stop).tab();
			if(r.known){buffer.append(r.fixed);}
			buffer.tab().append(r.partial).tab().append(r.known ? 1 : 0).tab().append(r.trim ? 1 : 0).tab()
				.append((writer.invocationsOnly ? r.invoked : r.eligible) ? 1 : 0).tab();
			if(r.known){buffer.append(Float.floatToIntBits(r.gc));}
			buffer.tab();
			if(r.known){buffer.append(Float.floatToIntBits(r.acceptedIdentity));}
			buffer.tab();
			if(r.known){buffer.append(Float.floatToIntBits(r.identity));}
			buffer.tab();
			if(r.known){buffer.append(r.quantum ? 1 : 0);}
			buffer.tab().append(r.min).tab();
			NcrnaTextIO.appendMaximumLength(buffer, r.max);
			buffer.tab().append(r.refinedStart).tab().append(r.refinedStop).tab();
			if(r.selected){buffer.append(r.finalStart);}
			buffer.tab();
			if(r.selected){buffer.append(r.finalStop);}
			buffer.tab();
			if(r.selected){buffer.append(r.finalIndex);}
			buffer.tab().append(r.family.positional ? "1-HOT" : "SLIDING").tab().append(r.family.radius).tab().append(r.scoredMask).tab()
				.append(r.returned ? 1 : 0).tab().append(r.selected ? 1 : 0).tab()
				.append(r.selected ? "SELECTED" : r.returned ? "PATH_REJECTED" : "SCAVENGER_REJECTED").tab();
			if(r.known){buffer.append(r.windowStart);}
			buffer.tab();
			if(r.known){buffer.append(r.windowStop);}
			buffer.tab();
			if(r.known){buffer.append(r.alignedLength);}
			buffer.nl();
			for(int i=0; i<r.count; i++){
				final Window w=r.windows[i];
				buffer.append("S\t").append(ordinal).tab().append(thread).tab().append(r.id).tab().append(w.end).tab()
					.append(w.context).tab().append(w.start).tab().append(w.stop).tab();
				for(int j=0; j<w.scores.length; j++){
					if(j>0){buffer.append(',');}
					if(w.available==null || w.available[j]){buffer.append(Float.floatToIntBits(w.scores[j]));}
				}
				buffer.nl();
			}
			synchronized(writer.output){writer.output.print(buffer);}
		}
		private final Writer writer;
		private final int ordinal, thread;
		private final IdentityHashMap<Orf, Record> pending=new IdentityHashMap<Orf, Record>();
		private final ArrayList<Record> records=new ArrayList<Record>();
		private final ArrayList<Family> families=new ArrayList<Family>();
		private final ByteBuilder buffer=new ByteBuilder();
		private String fullName, contig;
		private int length;
		private long nextID, readID;
	}

	/** Per-family worker adapter; active state cannot be shared across caller workers. */
	public static final class Family implements NcrnaBoundaryScorer.ScoreSink{
		private Family(final Worker owner_, final NcrnaFamily family_){
			owner=owner_;
			family=family_;
			positional=family.boundaryRefinerFactory instanceof NcrnaPositionalBoundaryRefiner.Resources;
			radius=owner.writer.invocationsOnly ? 0 : positional ? ((NcrnaPositionalBoundaryRefiner.Resources)family.boundaryRefinerFactory).radius()
				: NcrnaBoundaryScorer.centeredRadius(family.boundaryStartOffsets);
		}
		public void begin(final Orf orf, final int model, final boolean trim, final boolean eligible, final int fixed,
				final float gc, final float acceptedIdentity, final float identity, final boolean quantum, final int min, final int max,
				final int windowStart, final int windowStop, final int alignedLength){
			require(current==null && Float.isFinite(gc) && Float.isFinite(acceptedIdentity) && Float.isFinite(identity)
				&& min>0 && max>=min && (fixed&~3)==0,
				"Trace boundary seams cannot nest or synthesize feature/length metadata");
			require(windowStart>=0 && windowStop>=windowStart && windowStop<owner.length && alignedLength>0,
				"Eligibility provenance must retain the actual accepted window and aligned length");
			current=owner.add(this, orf, model, true, trim, eligible, fixed, gc, acceptedIdentity, identity,
				quantum, min, max, windowStart, windowStop, alignedLength);
		}
		public void complete(final Orf orf){
			complete(orf, orf.start, orf.stop);
		}
		/** Saves the effective report even when selection bounds have already been restored. */
		public void complete(final Orf orf, final int reportedStart, final int reportedStop){
			require(current!=null && owner.pending.get(orf)==current, "Complete the exact candidate whose boundary seam began");
			require(reportedStart>=0 && reportedStop>=reportedStart && reportedStop<owner.length,
				"Effective report must retain physical inclusive coordinates before strand conversion");
			current.refinedStart=reportedStart;
			current.refinedStop=reportedStop;
			current.complete=true;
			current=null;
			windowOffset=0;
		}
		/** Records the actual post-trim refiner call, rather than merely its eligibility predicate. */
		public void invoked(final Orf orf){
			require(invocationsOnly() && current!=null && owner.pending.get(orf)==current
				&& current.eligible && !current.invoked,
				"Invocation evidence must identify one eligible candidate immediately before its refiner call");
			current.invoked=true;
		}
		public boolean invocationsOnly(){return owner.writer.invocationsOnly;}
		/** Reuses the existing scavenger-return event, retaining unknown eligibility explicitly. */
		public void returned(final Orf orf){
			Record record=owner.pending.get(orf);
			if(record==null){
				int model=-1;
				for(int i=0; i<family.modelNames.length; i++){
					if(family.modelNames[i].equals(orf.trnaModel)){
						model=i;
						break;
					}
				}
				record=owner.add(this, orf, model, false, false, false, 0, 0, 0, 0,
					false, family.minLen, family.maxLen, -1, -1, -1);
				record.refinedStart=orf.start;
				record.refinedStop=orf.stop;
				record.complete=true;
			}
			require(record.complete && !record.returned, "Each traced object may return from its scavenger only once");
			record.returned=true;
		}
		/** Converts scalar window-relative centers back to original sense coordinates. */
		public void windowOffset(final int offset){windowOffset=offset;}
		@Override
		public void scores(final int end, final int context, final int start, final int stop,
				final float[] values, final boolean[] available){
			if(invocationsOnly()){return;}
			require(current!=null && current.eligible && end>=0 && end<2 && current.count<current.windows.length
				&& values!=null && values.length==2*radius+1, "Trace scores must belong to a real eligible centered endpoint window");
			require(available==null || available.length==values.length && available[radius],
				"Trace availability must cover the exact score window and retain its unchanged position");
			for(int i=0; i<values.length; i++){
				require(available!=null && !available[i] || Float.isFinite(values[i]),
					"A computed trace score must be finite before capture: index="+i);
			}
			current.windows[current.count++]=new Window(end, context, start+windowOffset, stop+windowOffset, values, available);
			current.scoredMask|=1<<end;
		}
		private final Worker owner;
		private final NcrnaFamily family;
		private final boolean positional;
		private final int radius;
		private Record current;
		private int windowOffset;
	}

	/** Owns the transferred numeric output; network buffers were copied by their producer. */
	private static final class Window{
		Window(final int end_, final int context_, final int start_, final int stop_, final float[] scores_, final boolean[] available_){
			end=end_;
			context=context_;
			start=start_;
			stop=stop_;
			scores=scores_;
			available=available_;
		}
		final int end, context, start, stop;
		final float[] scores;
		final boolean[] available;
	}

	private static final class Record{
		Record(final Family family_, final Orf orf, final int model_, final long id_, final boolean known_,
				final boolean trim_, final boolean eligible_, final int fixed_, final float gc_, final float acceptedIdentity_, final float identity_,
				final boolean quantum_, final int min_, final int max_, final int windowStart_, final int windowStop_, final int alignedLength_){
			family=family_;
			id=id_;
			model=NcrnaSupport.recordName(family.family.modelNames[model_]);
			strand=orf.strand;
			start=orf.start;
			stop=orf.stop;
			partial=orf.ncrnaPartial;
			known=known_;
			trim=trim_;
			eligible=eligible_;
			fixed=fixed_;
			gc=gc_;
			acceptedIdentity=acceptedIdentity_;
			identity=identity_;
			quantum=quantum_;
			min=min_;
			max=max_;
			windowStart=windowStart_;
			windowStop=windowStop_;
			alignedLength=alignedLength_;
		}
		long endpoints(final float margin){
			int five=start, three=stop;
			if(eligible && count>0){
				if(family.positional){
					for(int i=0; i<count; i++){
						final Window w=windows[i];
						final int offset=selected(w.scores, w.available, family.radius, margin, true);
						if(w.end==0){five+=offset;}else{three+=offset;}
					}
					final long span=(long)three-five+1;
					if(five<0 || three>=family.owner.length || span<min || span>max){
						five=start;
						three=stop;
					}
				}else{
					final Window first=windows[0];
					require(first.context==FIRST, "Scalar replay needs its actual first-end order");
					final int offset=selected(first.scores, first.available, family.radius, margin, false);
					if(first.end==0){five+=offset;}else{three+=offset;}
					boolean secondFound=false;
					for(int i=1; i<count; i++){
						final Window w=windows[i];
						if(w.context==(offset==0 ? HELD : MOVED)){
							final int second=selected(w.scores, w.available, family.radius, margin, false);
							if(w.end==0){five+=second;}else{three+=second;}
							secondFound=true;
							break;
						}
					}
					require(fixed!=0 || secondFound, "A two-end scalar replay must retain the needed held/moved second-end context");
				}
			}
			return ((long)five<<32)|(three&0xffffffffL);
		}
		final Family family;
		final long id;
		final String model;
		final int strand, start, stop, partial, fixed, min, max, windowStart, windowStop, alignedLength;
		final boolean known, trim, eligible, quantum;
		final float gc, acceptedIdentity, identity;
		final Window[] windows=new Window[3];
		int count, scoredMask, refinedStart, refinedStop, finalStart=-1, finalStop=-1, finalIndex=-1;
		boolean complete, returned, selected, invoked;
	}

	private static int selected(final float[] scores, final boolean[] available, final int radius, final float margin, final boolean positional){
		int best=-1;
		for(int i=0; i<scores.length; i++){
			if(available!=null && !available[i]){continue;}
			require(Float.isFinite(scores[i]), "Nonfinite computed trace scores cannot establish an endpoint choice");
			if(best<0 || scores[i]>scores[best]){best=i;}
		}
		require(best>=0 && Float.isFinite(scores[radius]), "A scored trace window must retain its actual current position");
		final int legacy=positional || scores[best]>scores[radius] ? best : radius;
		return NcrnaEndpointConfidence.select(scores, radius, legacy, margin, available)-radius;
	}
	private static boolean matches(final String family, final String target){
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family);
		return family.equals(target) || entry.group!=null && entry.group.id.equals(target);
	}
	/** Checks the explicit target before opening its output writer. */
	public static void validateTarget(final String target, final ArrayList<NcrnaFamily> families){
		validateTarget(target, families, false);
	}
	public static void validateTarget(final String target, final ArrayList<NcrnaFamily> families, final boolean invocationsOnly){
		boolean found=false;
		for(final NcrnaFamily family:families){
			if(matches(family.name, target)){
				validateFamily(family, invocationsOnly);
				found=true;
			}
		}
		require(found, "Boundary trace target is not enabled: "+target);
	}
	private static void validateFamily(final NcrnaFamily family, final boolean invocationsOnly){
		//TODO: Probable bug - raw/deferred glycine reporting can reject an NN proposal
		// after this capture, violating finish's transformed-refinement check. Reject
		// boundaryOnRawEndpoints or define that reporting contract before supporting
		// this separate route; the current alignment-centred cobalamin/SAM path is tested.
		NcrnaFamilyRegistry.sharedBoundaryMember(family.name);
		final boolean positional=family.boundaryRefinerFactory instanceof NcrnaPositionalBoundaryRefiner.Resources;
		final boolean scalar=family.boundary5NetTemplate!=null && family.boundary3NetTemplate!=null;
		if(invocationsOnly){
			require((positional || family.boundaryRefinerFactory==null && !scalar) && !family.voteEnds
				&& !family.nnCentre5Vote && !family.nnCentre3Vote && family.cdsStopConstraintFactory==null
				&& family.rrnaEndpointFeatures==null && family.boundaryNetsByModel==null,
				"Invocation-only evidence supports ordinary positional refiners or an unchanged member");
			return;
		}
		require(positional!=scalar && family.cdsStopConstraintFactory==null && family.rrnaEndpointFeatures==null
			&& (family.boundaryRefinerFactory==null || positional)
			&& !family.voteEnds && !family.nnCentre5Vote && !family.nnCentre3Vote
			&& !family.nnCutoffEnabled && Float.isFinite(family.nnCutoff),
			"Trace only ordinary ungated scalar or positional inference, without votes or another endpoint route");
		if(positional){
			require(!((NcrnaPositionalBoundaryRefiner.Resources)family.boundaryRefinerFactory).fullErrors(),
				"Boundary trace v1 does not support full-error veto heads or strict zero-margin selection");
			require(((NcrnaPositionalBoundaryRefiner.Resources)family.boundaryRefinerFactory).margin()==0,
				"Capture scores in the already-required ungated positional run");
		}else{
			final int radius=NcrnaBoundaryScorer.centeredRadius(family.boundaryStartOffsets);
			require(radius>0 && NcrnaBoundaryScorer.centeredRadius(family.boundaryStopOffsets)==radius
				&& family.boundaryMarginStart==0 && family.boundaryMarginStop==0 && family.boundaryNetsByModel==null,
				"Scalar trace needs the declared common centered radius and zero margins");
		}
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalStateException(why);}}

	public static final int FIRST=0, HELD=1, MOVED=2, INDEPENDENT=3;
	public static final String HEADER="#ncrna_boundary_trace\t3\n"
		+"#Coordinates: zero-based inclusive; anchors/windows/refined are original sense frame, final is genomic.\n"
		+"#Contexts: 0=first scalar end, 1=second with first held, 2=second with first moved, 3=independent positional end.\n"
		+"#Scores: index i is correction i-radius; present float fields use Float.intBitsToFloat; unavailable positions have empty CSV cells.\n"
		+"#Optional C fields: fixed/features/window metadata are empty when known=0; final bounds/index are empty when selected=0.\n"
		+"#Length: max_length is empty when the family has no upper length bound.\n"
		+"#I\tinput_ordinal0\tassembly\n"
		+"#C\tinput_ordinal0\tworker\ttrace_id\tread_id0\tcontig\tstrand\tfamily\tmodel\tcontig_length\tanchor_start0\tanchor_stop0\tfixed_mask\tpartial_mask\tknown\ttrim_succeeded\tnn_invoked\tgc_bits\taccepted_identity_bits\tidentity_bits\tquantum\tmin_length\tmax_length\trefined_start0\trefined_stop0\tfinal_start0\tfinal_stop0\tfinal_index0\tkind\tradius\tscored_end_mask\tscavenger_returned\tselected\toutcome\taccepted_window_start0\taccepted_window_stop0\taccepted_aligned_length\n"
		+"#S\tinput_ordinal0\tworker\ttrace_id\tend\tcontext\tcenter_start0\tcenter_stop0\tscore_bits_csv\n"
		+"#T\tinput_ordinal0\tworker\tread_id0\tcontig\tlength\tcandidates\tscavenger_returned\tselected\tunknown\teligible_without_scores\n";
	/** Same source/candidate columns, but nn_invoked means an observed custom-refiner call.
	 * No S rows or replay ability; zero scored_end_mask is expected, not missing inference. */
	public static final String INVOCATION_HEADER=HEADER.replace("#ncrna_boundary_trace\t3\n",
		"#ncrna_boundary_invocations\t3\n#Invocation-only: nn_invoked observes the custom refiner call; no score replay.\n");
	private NcrnaBoundaryTrace(){}
}
