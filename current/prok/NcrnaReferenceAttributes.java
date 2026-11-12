package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaHistoricalEvidence.Attributes;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/** Joins current original-GFF attributes onto all saved retained source-row keys.
 * Matches the saved target/type/coordinates/strand/ID and Parent/RFAM claims.
 * Other attributes are captured now, not certified as identical to uncaptured
 * historical bytes. Each distinct GFF path is read once, without FASTA observations.
 * Workers own distinct files and reusable parsers/buffers; metadata is read-only
 * after loading. Native writer jobs transfer complete buffers. No gene admission.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceAttributes implements Accumulator<NcrnaReferenceAttributes.Worker>{

	static void process(final String input, final String output, final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceProfiles profiles=NcrnaReferenceProfiles.load(input, catalog, roster, true);
		final NcrnaReferenceAttributes join=new NcrnaReferenceAttributes(catalog, profiles);
		join.loadHistory(input+".historical.tsv.gz");
		join.loadRelations(input+".relations.tsv.gz");
		try(Output out=new Output(output+".attributes.tsv.gz")){
			join.writer=out.writer;
			join.writer.addJob(new ByteBuilder().append("input_id\tgff\tline\traw_attributes\n"));
			join.runWorkers();
		}
		join.writeSummary(output);
	}

	private NcrnaReferenceAttributes(final NcrnaEvidenceCatalog catalog_, final NcrnaReferenceProfiles profiles_){
		catalog=catalog_;
		profiles=profiles_;
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			inputs.put(genome.id, new InputInfo(genome, i, profiles.selected(genome.id)));
			if(!profiles.selected(genome.id)){continue;}
			for(int j=0; j<genome.memberships(); j++){
				final String path=genome.membership(j).gff;
				final byte[] encoded=BuildNcrnaReference.field(b.clear(), path).toBytes();
				final String key=new String(encoded, StandardCharsets.UTF_8);
				if(!sources.contains(key)){sources.put(key, new SourceFile(path, encoded));}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved bindings       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadHistory(final String path){
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		InputInfo input=null;
		SourceFile source=null;
		long total=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, profiles.batchEvidence() ? BuildNcrnaReference.BATCH_HISTORICAL_HEADER : BuildNcrnaReference.HISTORICAL_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==(profiles.batchEvidence() ? 20 : 19), "Saved historical schema differs from its evidence generation");
				final InputInfo priorInput=input;
				final SourceFile priorSource=source;
				if(input==null || !p.termEquals(input.genome.id, 0)){input=inputs.get(p.parseString(0));}
				if(source==null || !same(p, 1, source.encoded)){source=sources.get(text(p, 1));}
				require(input!=null && source!=null && input.selected, "Historical row needs a selected catalog input and original GFF claim");
				if(input!=priorInput || source!=priorSource){
					boolean claimed=false;
					for(int i=0; i<input.genome.memberships(); i++){claimed|=source.path.equals(input.genome.membership(i).gff);}
					require(claimed && new File(source.path).isAbsolute(), "Historical source is not an absolute GFF claimed by this input");
				}
				final long line=p.parseLong(2);
				final int parents=p.parseInt(13), rfam=p.parseInt(14);
				require(line>0 && parents>=0 && rfam>=0, "Invalid saved source line or relationship counts");
				b.clear();
				for(int field=3; field<=7; field++){p.appendTerm(b, field).tab();}
				final int idOffset=b.length();
				p.appendTerm(b, 8);
				if(source.rows.isEmpty()){files.add(source);}
				source.rows.add(new Expected(input, line, b.toBytes(), idOffset, parents, rfam));
				input.rows++;
				total++;
			}
		}
		require(total==profiles.retained(), "Historical rows must conserve the saved retained-source total");
		for(int i=0; i<catalog.size(); i++){
			final InputInfo inputInfo=inputs.get(catalog.get(i).id);
			require(inputInfo.rows==profiles.retained(inputInfo.genome.id), "Historical row counts disagree for input "+inputInfo.genome.id);
		}
		for(final SourceFile sourceFile:files){
			Collections.sort(sourceFile.rows, ROW_ORDER);
			Expected previous=null;
			for(final Expected expected:sourceFile.rows){
				require(previous==null || ROW_ORDER.compare(previous, expected)!=0, "Duplicate input/GFF/line in saved history: "+sourceFile.path);
				previous=expected;
			}
		}
		Collections.sort(files, (a, z)->a.path.compareTo(z.path));
	}

	private void loadRelations(final String path){
		final LineParser1 p=new LineParser1('\t');
		InputInfo input=null;
		SourceFile source=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.RELATION_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Saved relationships require seven fields");
				// Cycle witnesses can describe unselected source rows; they remain in the original table.
				if(p.termEquals("CYCLE", 3)){continue;}
				final boolean parent=p.termEquals("PARENT", 3);
				require(parent || p.termEquals("RFAM", 3), "Unknown saved relationship kind");
				if(input==null || !p.termEquals(input.genome.id, 0)){input=inputs.get(p.parseString(0));}
				if(source==null || !same(p, 1, source.encoded)){source=sources.get(text(p, 1));}
				require(input!=null && source!=null, "Relationship must identify a catalog input and source");
				final Expected expected=source.find(p.parseLong(2), input.index);
				require(expected!=null, "Parent/RFAM relationship has no retained source row");
				p.setBounds(4);
				require(same(p.line(), p.a(), p.b(), expected.identity, expected.idOffset, expected.identity.length), "Relationship subject ID differs from its saved historical row");
				if(parent){
					require(expected.parentsSeen<expected.parents.length, "Too many saved Parent claims for a source row");
					expected.parents[expected.parentsSeen++]=p.parseByteArray(5);
				}else{
					require(expected.rfamSeen<expected.rfam.length, "Too many saved RFAM claims for a source row");
					expected.rfam[expected.rfamSeen++]=p.parseByteArray(5);
				}
			}
		}
		for(final SourceFile file:files){
			for(final Expected expected:file.rows){
				require(expected.parentsSeen==expected.parents.length && expected.rfamSeen==expected.rfam.length,
					"Missing saved Parent/RFAM claims at "+file.path+":"+expected.line);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker lifecycle     ----------------*/
	/*--------------------------------------------------------------*/

	private void runWorkers(){
		if(files.isEmpty()){return;}
		final int count=Math.min(Shared.threads(), files.size());
		final ArrayList<Worker> workers=new ArrayList<Worker>(count);
		for(int i=0; i<count; i++){workers.add(new Worker(this, i, count));}
		try{ThreadWaiter.startAndWait(workers, this);}
		catch(RuntimeException | Error failure){
			stop=true;
			final ArrayList<Worker> started=new ArrayList<Worker>();
			for(final Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
			ThreadWaiter.waitForThreadsToFinish(started);
			throw failure;
		}
		for(final Worker worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Original-GFF attribute join failed", worker.failure);}
		}
		require(!stop && completed==files.size() && emitted==profiles.retained(), "All original source bindings must match before completion");
	}
	@Override
	public void accumulate(final Worker worker){
		assert(worker.getState()==Thread.State.TERMINATED) : "Only joined workers surrender their sole-owned counters";
		completed+=worker.completed;
		emitted+=worker.emitted;
	}
	@Override
	public boolean success(){return !stop;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	static final class Worker extends Thread{
		Worker(final NcrnaReferenceAttributes owner_, final int first_, final int stride_){
			owner=owner_;
			first=first_;
			stride=stride_;
		}
		@Override
		public void run(){
			try{
				for(int i=first; i<owner.files.size() && !owner.stop; i+=stride){read(owner.files.get(i));}
				flush();
			}catch(Throwable problem){
				failure=problem;
				owner.stop=true;
			}
		}
		private void read(final SourceFile file){
			int next=0;
			byte[] previousTarget=null;
			String target=null;
			try(Input in=new Input(file.path)){
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					if(owner.stop){return;}
					file.lines++;
					if(EukOrganellePanelInventory.starts(row, "##FASTA") && (row.length==7 || row[7]<=32)){break;}
					if(next==file.rows.size() || file.rows.get(next).line!=file.lines){continue;}
					p.set(row);
					require(p.terms()==9 && row.length>0 && row[0]!='#', "Saved line no longer identifies a GFF feature");
					if(previousTarget==null || !same(p, 0, previousTarget)){
						previousTarget=p.parseByteArray(0);
						target=attributes.decode(previousTarget, 0, previousTarget.length);
					}
					attributes.read(p);
					BuildNcrnaReference.field(identity.clear(), target).tab();
					BuildNcrnaReference.field(identity, NcrnaHistoricalEvidence.type(p)).tab();
					identity.append(EukPanelGffValidation.coordinate(p, 3)).tab().append(EukPanelGffValidation.coordinate(p, 4)).tab();
					p.setBounds(6);
					identity.append(p.b()==p.a()+1 ? (char)row[p.a()] : '?').tab();
					BuildNcrnaReference.field(identity, attributes.id==null ? "NA" : attributes.id);
					file.uniqueLines++;
					while(next<file.rows.size() && file.rows.get(next).line==file.lines){
						final Expected expected=file.rows.get(next++);
						require(same(identity.array, 0, identity.length(), expected.identity, 0, expected.identity.length),
							"Current target/type/coordinates/strand/ID differ from the saved source row");
						checkClaims(attributes.parents, expected.parents, "Parent");
						checkClaims(attributes.rfam, expected.rfam, "RFAM");
						buffer.append(expected.input.genome.id).tab().append(file.encoded).tab().append(expected.line).tab();
						p.setBounds(8);
						BuildNcrnaReference.field(buffer, row, p.a(), p.b()).nl();
						emitted++;
						file.emitted++;
						if(buffer.length()>=BLOCK_SIZE){flush();}
					}
				}
				require(next==file.rows.size(), "Original GFF ended before all saved retained lines were matched");
			}catch(RuntimeException problem){
				throw new IllegalArgumentException(file.path+":"+file.lines+": "+problem.getMessage(), problem);
			}
			completed++;
		}
		private void checkClaims(final String[] actual, final byte[][] expected, final String kind){
			if(actual.length!=expected.length){throw new IllegalArgumentException("Current "+kind+" count differs from saved relationships");}
			for(int i=0; i<actual.length; i++){
				BuildNcrnaReference.field(claim.clear(), actual[i]);
				if(!same(claim.array, 0, claim.length(), expected[i], 0, expected[i].length)){
					throw new IllegalArgumentException("Current "+kind+" identity/order differs from saved relationships");
				}
			}
		}
		private void flush(){
			if(buffer.length()==0){return;}
			owner.writer.addJob(buffer);
			buffer=new ByteBuilder(BLOCK_SIZE);
		}
		final NcrnaReferenceAttributes owner;
		final int first, stride;
		final LineParser1 p=new LineParser1('\t');
		final Attributes attributes=new Attributes();
		final ByteBuilder identity=new ByteBuilder(), claim=new ByteBuilder();
		ByteBuilder buffer=new ByteBuilder(BLOCK_SIZE);
		long emitted;
		int completed;
		Throwable failure;
	}

	private void writeSummary(final String prefix){
		final ByteBuilder b=new ByteBuilder();
		long lines=0, selected=0;
		try(Output out=new Output(prefix+".attribute_sources.tsv")){
			out.writer.println("gff\tsaved_input_rows\tmatched_input_rows\tunique_selected_lines\tsource_lines_read");
			for(final SourceFile file:files){
				out.writer.print(b.clear().append(file.encoded).tab().append(file.rows.size()).tab().append(file.emitted).tab().append(file.uniqueLines).tab().append(file.lines).nl());
				lines+=file.lines;
				selected+=file.uniqueLines;
			}
		}
		try(Output out=new Output(prefix+".attribute_summary.tsv")){
			out.writer.println("inputs_observed\tcatalog_inputs\tretained_source_rows\tmatched_attribute_rows\tsource_paths_read\tunique_selected_lines\tsource_lines_read\tevidence_status\tstatus");
			b.clear().append(profiles.observed()).tab().append(catalog.size()).tab().append(profiles.retained()).tab().append(emitted).tab();
			out.writer.print(b.append(files.size()).tab().append(selected).tab().append(lines).tab().append(profiles.evidenceState()).append("\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class InputInfo{
		InputInfo(final Genome genome_, final int index_, final boolean selected_){
			genome=genome_;
			index=index_;
			selected=selected_;
		}
		final Genome genome;
		final int index;
		final boolean selected;
		long rows;
	}
	private static final class SourceFile{
		SourceFile(final String path_, final byte[] encoded_){
			path=path_;
			encoded=encoded_;
		}
		Expected find(final long line, final int input){
			int low=0, high=rows.size()-1;
			while(low<=high){
				final int mid=(low+high)>>>1;
				final Expected row=rows.get(mid);
				final int c=row.line==line ? Integer.compare(row.input.index, input) : Long.compare(row.line, line);
				if(c<0){low=mid+1;}
				else if(c>0){high=mid-1;}
				else{return row;}
			}
			return null;
		}
		final String path;
		final byte[] encoded;
		final ArrayList<Expected> rows=new ArrayList<Expected>();
		long lines, uniqueLines, emitted;
	}
	private static final class Expected{
		Expected(final InputInfo input_, final long line_, final byte[] identity_, final int idOffset_, final int parentCount, final int rfamCount){
			input=input_;
			line=line_;
			identity=identity_;
			idOffset=idOffset_;
			parents=parentCount==0 ? EMPTY : new byte[parentCount][];
			rfam=rfamCount==0 ? EMPTY : new byte[rfamCount][];
		}
		final InputInfo input;
		final long line;
		final byte[] identity;
		final int idOffset;
		final byte[][] parents, rfam;
		int parentsSeen, rfamSeen;
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		return same(p.line(), p.a(), p.b(), expected, 0, expected.length);
	}
	private static boolean same(final byte[] a, final int from, final int to, final byte[] b, final int start, final int end){
		if(to-from!=end-start){return false;}
		for(int i=0; i<to-from; i++){if(a[from+i]!=b[start+i]){return false;}}
		return true;
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final Comparator<Expected> ROW_ORDER=(a, b)->a.line==b.line ? Integer.compare(a.input.index, b.input.index) : Long.compare(a.line, b.line);
	private static final byte[][] EMPTY=new byte[0][];
	static final String[] OUTPUT_SUFFIXES={".attributes.tsv.gz", ".attribute_sources.tsv", ".attribute_summary.tsv"};
	private static final int BLOCK_SIZE=65536;
	private final NcrnaEvidenceCatalog catalog;
	private final NcrnaReferenceProfiles profiles;
	private final ObjectMap<String,InputInfo> inputs=new ObjectMap<String,InputInfo>(String.class, InputInfo.class);
	private final ObjectMap<String,SourceFile> sources=new ObjectMap<String,SourceFile>(String.class, SourceFile.class);
	private final ArrayList<SourceFile> files=new ArrayList<SourceFile>();
	private ByteStreamWriter writer;
	private volatile boolean stop;
	private long emitted;
	private int completed;
}
