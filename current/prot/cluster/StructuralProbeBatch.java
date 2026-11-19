package prot.cluster;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prot.ProteinSequence;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Persisted independent probe requests and compact native results.
 * One engine/context serves a batch, but every row edits the same original population.
 * The caller preflights every output path before invoking run; no scientific sampling
 * or objective-based candidate filtering occurs here.
 * @author Keqing
 */
public final class StructuralProbeBatch{

	/** Reads explicit candidate identities. A multi-row batch needs # in its codec output pattern. */
	static Plan load(final String path,final String outputPattern){
		if(path==null || outputPattern==null){throw new IllegalArgumentException("Probe batches need a request table and codec output pattern");}
		return read(path,outputPattern);
	}

	/** Local evaluation consumes the same persisted requests without exporting model overlays. */
	static Plan loadRequests(final String path){
		if(path==null){throw new IllegalArgumentException("Local probes need a persisted request table");}
		return read(path,null);
	}

	private static Plan read(final String path,final String outputPattern){
		final ArrayList<Request> rows=new ArrayList<Request>();
		final IntHashMap ids=new IntHashMap();
		final ByteFile reader=ByteFile.makeByteFile(path,false);
		try{
			final byte[] header=reader.nextLine();
			if(header==null || !REQUEST_HEADER.equals(new String(header))){throw new IllegalArgumentException("Unexpected probe request header");}
			final LineParser1 fields=new LineParser1('\t');
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				fields.set(line);
				if(fields.terms()!=5){throw new IllegalArgumentException("Probe requests need exactly five columns");}
				final int id=identifier(fields,0),left=identifier(fields,2),right=identifier(fields,3);
				final LayerCakeSchedule.Phase phase;
				if(fields.termEquals("fission",1)){phase=LayerCakeSchedule.Phase.FISSION;}
				else if(fields.termEquals("fusion",1)){phase=LayerCakeSchedule.Phase.FUSION;}
				else{throw new IllegalArgumentException("Probe operation must be fission or fusion");}
				if(id<0 || ids.get(id)>=0 || left<0 || (phase==LayerCakeSchedule.Phase.FISSION ? right!=-1 : right<0 || right==left)){
					throw new IllegalArgumentException("Probe IDs must be unique/nonnegative and endpoints must describe one parent or distinct pair");
				}
				ids.put(id,1);
				final String out=outputPattern==null ? null : outputPattern.replace("#",Integer.toString(id));
				if(out!=null && (out.indexOf('\t')>=0 || out.indexOf('\n')>=0 || out.indexOf('\r')>=0)){
					throw new IllegalArgumentException("Probe codec paths must fit an unescaped TSV field");
				}
				rows.add(new Request(id,decimal(fields,4),phase,left,right,out));
			}
		}finally{if(reader.close()){throw new IllegalStateException("Failed probe request input: "+path);}}
		if(rows.isEmpty() || outputPattern!=null && rows.size()>1 && outputPattern.indexOf('#')<0){
			throw new IllegalArgumentException("Require a nonempty probe table and # in a multi-row codec output pattern");
		}
		return new Plan(rows.toArray(new Request[rows.size()]));
	}

	/** Request metadata must not inherit Parse.parseLong's unchecked wraparound. */
	private static long decimal(final LineParser1 fields,final int term){
		fields.setBounds(term);
		final byte[] line=fields.line();int at=fields.a();final int end=fields.b();
		if(at>=end){throw new IllegalArgumentException("Empty probe numeric field");}
		final boolean negative=line[at]=='-';if(negative){at++;}
		if(at>=end){throw new IllegalArgumentException("Missing probe numeric digits");}
		final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
		long value=0;
		for(; at<end; at++){
			final int digit=line[at]-'0';
			if(digit<0 || digit>9 || value<limit/10){throw new IllegalArgumentException("Invalid or overflowing probe numeric field");}
			value*=10;
			if(value<limit+digit){throw new IllegalArgumentException("Overflowing probe numeric field");}
			value-=digit;
		}
		return negative ? value : -value;
	}
	private static int identifier(final LineParser1 fields,final int term){
		final long value=decimal(fields,term);
		if(value< -1 || value>=Integer.MAX_VALUE){throw new IllegalArgumentException("Probe identifier is outside the supported range");}
		return (int)value;
	}

	/** Runs only after the common original-member loader and complete path preflight.
	 * Rejected rows write an empty codec and no lineage, so overwrite mode cannot leave
	 * an older successful payload masquerading as a new result. The consumer rejects
	 * that empty overlay; callers must also require feasible=1 from the result table.
	 */
	static void run(final Plan plan,final LayerCakeClustering engine,final Family[] originals,
			final ProteinSequence[] sequences,final ClusterSketchCodec.Header header,final double temperature,
			final String report,final String workReport,final String members,final boolean overwrite){
		plan.validate(originals.length);
		final ByteStreamWriter[] writers=new ByteStreamWriter[3];
		try{
			writers[0]=open(report,overwrite);writers[1]=open(workReport,overwrite);writers[2]=open(members,overwrite);
			writers[0].println(RESULT_HEADER);writers[1].println(WORK_HEADER);
			writers[2].println("candidate_index\tparent_input_index\tinput_id\tfamily_id");
			final ByteBuilder row=new ByteBuilder();final long[] before=new long[8];
			for(Request request:plan.rows){
				counters(engine,before);
				final long start=System.nanoTime();
				final LayerCakeClustering.Probe result=engine.probe(request.phase,request.left,request.right,temperature);
				final long elapsed=System.nanoTime()-start;
				final Family[] replacements=result.feasible() ? result.replacements() : EMPTY;
				ClusterSketchCodec.write(request.output,overwrite,replacements,header.withFamilies(replacements.length));
				final long affected=(long)originals[result.inputA].size()+(result.inputB<0 ? 0 : originals[result.inputB].size());
				resultRow(row,request,result,affected);writers[0].println(row);
				row.clear().append(request.id).tab().append(elapsed).tab()
					.append(engine.qualityAlignments-before[0]).tab().append(engine.neighborAlignments-before[1]).tab()
					.append(engine.rebuildAlignments-before[2]).tab().append(engine.selectionAlignments-before[3]).tab()
					.append(engine.memberPostings-before[4]).tab().append(engine.neighborPostings-before[5]).tab()
					.append(engine.fissionRetainedEdges-before[6]).tab().append(engine.fissionUnionEdges-before[7]).tab().append("NA");
				writers[1].println(row);
				long written=0;
				for(Family family:replacements){
					for(int i=0; i<family.size(); i++){
						final int ordinal=family.memberIndex(i);
						row.clear().append(request.id).tab().append(ordinal).tab().append(sequences[ordinal].id).tab().append(family.familyId);
						writers[2].println(row);written++;
					}
				}
				if(result.feasible() && written!=affected){throw new IllegalStateException("Probe export must retain every affected original member");}
			}
		}finally{close(writers);}
	}

	/** Parent/union are the combined model; zero counts and -1 IDs denote the absent shape only. */
	private static void resultRow(final ByteBuilder row,final Request request,final LayerCakeClustering.Probe p,final long affected){
		row.clear().append(request.id).tab().append(request.seed).tab().append(request.phase.name().toLowerCase(Locale.ROOT))
			.tab().append(p.feasible() ? 1 : 0).tab().append(p.reason).tab().append(p.inputA).tab().append(p.inputB).tab().append(p.inputA);
		if(p.inputB>=0){row.append(',').append(p.inputB);}
		row.tab().append(request.output);
		if(p.feasible()){
			final boolean split=request.phase==LayerCakeSchedule.Phase.FISSION;
			row.tab().append(split ? p.parent.familyId : -1).tab().append(p.left.familyId).tab().append(p.right.familyId)
				.tab().append(split ? -1 : p.parent.familyId).tab().append(split ? p.parent.size() : 0)
				.tab().append(p.left.size()).tab().append(p.right.size()).tab().append(split ? 0 : p.parent.size());
		}else{missing(row,8);}
		row.tab().append(p.slotCost).tab().appendSlow(p.temperature);
		if(p.feasible()){
			final LayerCakeClustering.ProbeFeatures f=p.features;
			row.tab().appendSlow(f.balance).tab().appendSlow(f.cohesionDelta).tab().appendSlow(f.separationDelta)
				.tab().appendSlow(f.sizePressure).tab().appendSlow(f.slotDelta)
				.tab().appendSlow(f.parent.cohesion).tab().appendSlow(f.parent.neighbor)
				.tab().appendSlow(f.left.cohesion).tab().appendSlow(f.left.neighbor)
				.tab().appendSlow(f.right.cohesion).tab().appendSlow(f.right.neighbor);
		}else{missing(row,11);}
		row.tab().append(affected);
	}
	private static void missing(final ByteBuilder row,final int count){for(int i=0; i<count; i++){row.tab().append("NA");}}
	private static void counters(final LayerCakeClustering engine,final long[] into){
		into[0]=engine.qualityAlignments;into[1]=engine.neighborAlignments;into[2]=engine.rebuildAlignments;into[3]=engine.selectionAlignments;
		into[4]=engine.memberPostings;into[5]=engine.neighborPostings;into[6]=engine.fissionRetainedEdges;into[7]=engine.fissionUnionEdges;
	}
	private static ByteStreamWriter open(final String path,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);writer.start();return writer;
	}
	private static void close(final ByteStreamWriter[] writers){
		RuntimeException failure=null;
		for(ByteStreamWriter writer:writers){
			if(writer==null){continue;}
			try{if(writer.poisonAndWait()){throw new IllegalStateException("Failed probe batch output");}}
			catch(RuntimeException e){if(failure==null){failure=e;}else{failure.addSuppressed(e);}}
		}
		if(failure!=null){throw failure;}
	}

	static final class Plan{
		Plan(final Request[] rows_){rows=rows_;}
		int size(){return rows.length;}
		String[] outputs(){
			final String[] result=new String[rows.length];
			for(int i=0; i<rows.length; i++){
				if(rows[i].output==null){throw new IllegalStateException("Local request plans have no model-overlay output paths");}
				result[i]=rows[i].output;
			}
			return result;
		}
		void validate(final int familyCount){
			for(Request row:rows){if(row.left>=familyCount || row.right>=familyCount){throw new IllegalArgumentException("Probe endpoint is absent from the initial global roster");}}
		}
		final Request[] rows;
	}
	static final class Request{
		Request(final int id_,final long seed_,final LayerCakeSchedule.Phase phase_,final int left_,final int right_,final String out){
			id=id_;seed=seed_;phase=phase_;left=left_;right=right_;output=out;
		}
		final int id,left,right;final long seed;final LayerCakeSchedule.Phase phase;final String output;
	}
	public static final String REQUEST_HEADER="candidate_index\toperation\tleft_id\tright_id\tsample_seed";
	public static final String RESULT_HEADER="candidate_index\tsample_seed\toperation\tfeasible\treason\tinput_a\tinput_b\tremoval_ids\toverlay_path\tparent_id\tleft_id\tright_id\tunion_id\tparent_members\tleft_members\tright_members\tunion_members\tslot_cost\ttemperature\tbalance\tweighted_cohesion_delta\tweighted_separation_delta\tsize_pressure\tslot_delta\tcombined_quality\tcombined_neighbor\tleft_quality\tleft_neighbor\tright_quality\tright_neighbor\taffected_original_members";
	public static final String WORK_HEADER="candidate_index\tproposal_elapsed_ns\tquality_alignments\tneighbor_alignments\trebuild_alignments\tselection_alignments\tmember_postings\tneighbor_postings\tretained_edges\tunion_edges\tprocess_cpu_nanos";
	private static final Family[] EMPTY=new Family[0];
}
