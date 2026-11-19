package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Assembly-contribution census of the actual kept RF00174 construction members.
 * Uses the byte-file template's parser and native stream ownership. Reads saved
 * metadata only; this first stage does not assess exact gene-sequence overlap.
 * Known assembly aliases come exclusively from the accepted resolved group table.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlTrainingOverlap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlTrainingOverlap.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String genomes=null, groups=null, sample=null, census=null, captured=null;
			int inputCount=1000, genomeCount=39979, modelCount=508, memberCount=85741, orphanCount=15, smallCount=0;
			int capturedCount=1000;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("genomes")){genomes=b;}
				else if(a.equals("groups")){groups=b;}
				else if(a.equals("sample")){sample=b;}
				else if(a.equals("census")){census=b;}
				else if(a.equals("captured")){captured=b;}
				else if(a.equals("expectedcaptured")){capturedCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedinputs")){inputCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgenomes")){genomeCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmodels")){modelCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmembers")){memberCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedorphans")){orphanCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedsmall")){smallCount=Parse.parseIntKMG(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown training-overlap parameter: "+arg);}
			}
			require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use in=complete_assignments.tsv out=fresh_prefix without overwrite, append or sampling");
			process(parser.in1, genomes, groups, sample, census, parser.out1,
				inputCount, genomeCount, modelCount, memberCount, orphanCount, smallCount, captured, capturedCount);
			pp.outstream.println("NCRNA_DDL_ASSEMBLY_CONTRIBUTION_COMPLETE_EXACT_LOCUS_OVERLAP_NOT_ASSESSED");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Expected counts bind this census to the already verified voting population. */
	public static void process(final String assignments, final String genomes, final String groups,
			final String sample, final String census, final String output, final int inputCount,
			final int genomeCount, final int modelCount, final int memberCount, final int orphanCount, final int smallCount){
		process(assignments, genomes, groups, sample, census, output, inputCount, genomeCount, modelCount,
			memberCount, orphanCount, smallCount, null, 0);
	}

	/** Optional captured accounting selects one noncontributing representative per assembly.
	 * Already captured aliases take precedence; their native slots are preserved for reuse. */
	public static void process(final String assignments, final String genomes, final String groups,
			final String sample, final String census, final String output, final int inputCount,
			final int genomeCount, final int modelCount, final int memberCount, final int orphanCount, final int smallCount,
			final String captured, final int capturedCount){
		require(assignments!=null && genomes!=null && groups!=null && sample!=null && census!=null && output!=null
			&& inputCount>0 && genomeCount>0 && modelCount>0 && memberCount>0 && orphanCount>=0 && smallCount>=0,
			"Training overlap requires every saved metadata input and explicit complete positive population sizes");
		RrnaResourceIO.fresh(output, ".inputs.tsv", ".summary.tsv");
		if(captured!=null){RrnaResourceIO.fresh(output, ".reuse.tsv", ".uncaptured.selection.tsv", ".uncaptured.config", ".expansion.tsv");}
		final NcrnaDdlTrainingOverlap audit=new NcrnaDdlTrainingOverlap(inputCount, modelCount);
		audit.sample(sample);
		audit.groups(groups);
		if(captured!=null){audit.captured(captured, capturedCount);}
		audit.genomes(genomes, genomeCount);
		audit.census(census);
		audit.assignments(assignments, memberCount, orphanCount, smallCount);
		audit.write(output);
		if(captured!=null){audit.expansion(output);}
	}

	private NcrnaDdlTrainingOverlap(final int inputs, final int models){
		samples=new Sample[inputs];
		expected=new int[models];
		observed=new int[models];
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved Inputs         ----------------*/
	/*--------------------------------------------------------------*/

	private void sample(final String path){
		int count=0;
		try(Input in=new Input(path)){
			final byte[] header=in.reader.nextLine();
			selection=Arrays.equals(header, SELECTION_HEADER.getBytes(StandardCharsets.US_ASCII));
			require(selection || Arrays.equals(header, SAMPLE_HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Use the actual saved rank accounting or frozen eight-column benchmark selection schema");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(count<samples.length && (selection ? p.terms()==8 : p.terms()==13 && p.parseInt(0)==count),
					"Sample inputs must retain their complete actual selection or ordered saved accounting rows");
				final int id=p.parseInt(selection ? 0 : 1);
				final String fasta=text(selection ? 1 : 2);
				final long bases=p.parseLong(selection ? 2 : 4);
				final int records=p.parseInt(3);
				require(id>0 && fasta.startsWith("/") && fasta.indexOf(',')<0 && !sampleIds.contains(id) && samplePaths.add(fasta),
					"Sample IDs and literal paths must identify unique original evaluation inputs");
				require(bases>0 && records>0, "Frozen input sizes must describe real nonempty genomes, not placeholder accounting");
				samples[count]=new Sample(id, fasta, bases, records, selection ? row : null);
				sampleIds.put(id, count++);
			}
		}
		require(count==samples.length, "The assembly census must retain the complete declared sample, including zero-target inputs");
	}

	/** Bind the previous capture by original ID/path/size, never by its output order. */
	private void captured(final String path, final int expectedCount){
		require(selection && expectedCount>0, "Expansion needs the complete frozen selection and a declared prior captured population");
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, SAMPLE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==13 && p.parseInt(0)==count, "Prior capture accounting must retain its original dense slots");
				final int slot=sampleIds.get(p.parseInt(1));
				require(slot>=0, "Every captured input must occur in the complete evaluation selection");
				final Sample input=samples[slot];
				require(input.capturedSlot<0 && input.fasta.equals(text(2)) && input.records==p.parseInt(3) && input.bases==p.parseLong(4),
					"Prior capture ID, literal path and actual genome sizes must bind one frozen evaluation input each");
				input.capturedSlot=count++;
			}
		}
		require(count==expectedCount, "Expansion must account for the complete previous captured population");
	}

	private void groups(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==14 && (p.termEquals("evaluation", 0) || p.termEquals("discovery", 0)),
					"Assembly aliases must use the complete accepted cohort/group schema");
				final int group=p.parseInt(10);
				require(group>=0, "Resolved assembly group IDs must be nonnegative");
				final String fasta=text(2), canonical=text(3);
				putPath(fasta, group);
				putPath(canonical, group);
				require(p.termEquals("true", 12) || p.termEquals("false", 12), "Unresolved-pair status must be an explicit boolean");
				if(p.termEquals("true", 12)){uncertainGroups.put(group, 1);}
				if(p.termEquals("evaluation", 0)){
					final int slot=sampleIds.get(p.parseInt(1));
					if(slot<0){continue;}
					final Sample input=samples[slot];
					require(input.group<0 && input.fasta.equals(fasta),
						"Each sample input must bind one exact evaluation ID/path row, never a taxid guess");
					input.group=group;
				}
			}
		}
	}

	private void putPath(final String path, final int group){
		require(path.startsWith("/"), "Saved assembly binding requires absolute literal and canonical FASTA paths");
		final int old=pathGroups.get(path);
		require(old<0 || old==group, "A saved FASTA path cannot silently move between resolved assembly groups: "+path);
		pathGroups.put(path, group);
	}

	private void genomes(final String path, final int expectedCount){
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11 && p.parseInt(0)==count++, "The donor genome roster must preserve its complete ordinal order");
				final String genome=text(1), fasta=text(7);
				require(!genome.isEmpty() && !genomeGroups.contains(genome) && fasta.startsWith("/"),
					"Each donor source name needs one original genome-path binding");
				genomeGroups.put(genome, pathGroups.get(fasta));
			}
		}
		require(count==expectedCount, "Donor assembly binding must read every original genome, not only positives");
	}

	private void census(final String path){
		long members=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaMemberLengthBounds.HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==19 && p.termEquals("RF00174", 1) && modelIds.size()<expected.length,
					"The frozen census must retain every RF00174 model and its contribution count");
				final String model=text(0);
				final int count=p.parseInt(3), index=modelIds.size();
				require(count>0 && !modelIds.contains(model), "Frozen model IDs must be unique with nonempty member populations");
				modelIds.put(model, index);
				expected[index]=count;
				members+=count;
			}
		}
		require(modelIds.size()==expected.length && members>0 && members<=Integer.MAX_VALUE,
			"The complete frozen model census must define a nonempty bounded voting population");
		keptExpected=(int)members;
	}

	private void assignments(final String path, final int expectedCount, final int expectedOrphans, final int expectedSmall){
		final ObjectSet<String> headers=new ObjectSet<String>(String.class);
		int orphans=0, small=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, ASSIGNMENT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Every builder assignment must retain the complete record header and disposition");
				final String header=text(0);
				require(headers.add(header), "A construction record cannot contribute through duplicate assignment rows");
				if(p.termEquals("orphan", 5) || p.termEquals("small", 5)){
					require(p.termEquals("", 4), "Nonvoting records must not name a retained model");
					if(p.termEquals("orphan", 5)){orphans++;}else{small++;}
					continue;
				}
				require(p.termEquals("kept", 5), "Unknown builder dispositions must not acquire a vote");
				final int model=modelIds.get(RrnaResourceIO.first(text(4)));
				require(model>=0, "Every kept assignment must name one frozen current model");
				observed[model]++;
				kept++;
				final String genome=RrnaResourceIO.field(header, "source");
				require(genomeGroups.contains(genome), "A kept member's original source is absent from the complete genome roster: "+genome);
				donorGenomes.add(genome);
				final int group=genomeGroups.get(genome);
				if(group<0){unmappedVoters++; continue;}
				final int old=groupVoters.get(group);
				groupVoters.put(group, old<0 ? 1 : old+1);
				if(uncertainGroups.contains(group)){uncertainDonors.put(group, 1);}
			}
		}
		require(headers.size()==expectedCount && kept==keptExpected && orphans==expectedOrphans && small==expectedSmall
			&& (long)kept+orphans+small==expectedCount, "Kept and nonvoting assignments must conserve the frozen complete population");
		for(int i=0; i<expected.length; i++){
			require(observed[i]==expected[i], "Each model's actual kept assignment count must equal the frozen census; model="+i);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		final int[] inputs=new int[3];
		final IntHashMap[] groups={new IntHashMap(), new IntHashMap(), new IntHashMap()};
		int unmappedInputs=0;
		try(Output out=new Output(prefix+".inputs.tsv"); Output summary=new Output(prefix+".summary.tsv")){
			out.writer.println("input_slot\tinput_id\tfasta\tresolved_assembly_group\tkept_member_contributions\tassembly_status\tinterpretation");
			for(int i=0; i<samples.length; i++){
				final Sample sample=samples[i];
				final int votes=sample.group<0 ? -1 : groupVoters.get(sample.group);
				final int state=state(sample);
				inputs[state]++;
				if(sample.group<0){unmappedInputs++;}else{groups[state].put(sample.group, 1);}
				b.clear().append(i).tab().append(sample.id).tab().append(sample.fasta.getBytes(StandardCharsets.UTF_8))
					.tab().append(sample.group).tab().append(votes>0 ? votes : 0).tab().append(STATES[state])
					.append("\tASSEMBLY_CONTRIBUTION_ONLY_EXACT_SEQUENCE_AND_SOURCE_LOCUS_NOT_ASSESSED\n");
				out.writer.print(b);
			}
			require(inputs[0]+inputs[1]+inputs[2]==samples.length, "Seen, no-known-contribution and unknown input counts must conserve the sample");
			summary.writer.println("sample_inputs\tknown_contributor_inputs\tno_known_contributor_inputs\tunknown_inputs\tknown_contributor_groups\tno_known_contributor_groups\tunknown_resolved_groups\tunmapped_sample_inputs\tkept_members\tdonor_source_genomes\tunmapped_voters\tuncertain_donor_groups\tstatus");
			summary.writer.print(b.clear().append(samples.length).tab().append(inputs[0]).tab().append(inputs[1]).tab().append(inputs[2])
				.tab().append(groups[0].size()).tab().append(groups[1].size()).tab().append(groups[2].size()).tab().append(unmappedInputs)
				.tab().append(kept).tab().append(donorGenomes.size()).tab().append(unmappedVoters).tab().append(uncertainDonors.size())
				.append("\tKEPT_ASSIGNMENT_ASSEMBLY_CENSUS_NOT_EXACT_SEQUENCE_OR_HOLDOUT_PROOF\n"));
		}
	}

	private int state(final Sample sample){
		assert(sample!=null) : "Every contribution decision must retain its original sample binding";
		if(sample.group>=0 && groupVoters.get(sample.group)>0){return 0;}
		return sample.group<0 || unmappedVoters>0 || uncertainDonors.size()>0 || uncertainGroups.contains(sample.group) ? 2 : 1;
	}

	/** One representative per negative group; unknowns are never selected as negatives. */
	private void expansion(final String prefix){
		final IntHashMap representatives=new IntHashMap();
		for(int i=0; i<samples.length; i++){
			final Sample sample=samples[i];
			if(state(sample)!=1){continue;}
			final int previous=representatives.get(sample.group);
			if(previous<0 || sample.capturedSlot>=0 && (samples[previous].capturedSlot<0
				|| sample.capturedSlot<samples[previous].capturedSlot)){
				representatives.put(sample.group, i);
			}
		}
		int reused=0, additional=0;
		long bases=0, records=0;
		try(Output reuse=new Output(prefix+".reuse.tsv"); Output selected=new Output(prefix+".uncaptured.selection.tsv");
			Output config=new Output(prefix+".uncaptured.config"); Output summary=new Output(prefix+".expansion.tsv")){
			reuse.writer.println("resolved_assembly_group\tcaptured_slot\tinput_id\tfasta");
			selected.writer.println(SELECTION_HEADER);
			for(int i=0; i<samples.length; i++){
				final Sample sample=samples[i];
				if(sample.group<0 || representatives.get(sample.group)!=i){continue;}
				assert(state(sample)==1 && sample.selectionRow!=null)
					: "Expansion consumes actual frozen metadata only for resolved noncontributing assembly representatives";
				if(sample.capturedSlot>=0){
					reused++;
					reuse.writer.print(b.clear().append(sample.group).tab().append(sample.capturedSlot).tab().append(sample.id)
						.tab().append(sample.fasta.getBytes(StandardCharsets.UTF_8)).nl());
				}else{
					additional++;
					bases+=sample.bases;
					records+=sample.records;
					selected.writer.print(b.clear().append(sample.selectionRow).nl());
					config.writer.print(b.clear().append("in=").append(sample.fasta.getBytes(StandardCharsets.UTF_8)).nl());
				}
			}
			require(reused+additional==representatives.size(), "Reused and uncaptured groups must partition every eligible assembly exactly once");
			summary.writer.println("eligible_groups\treused_groups\tuncaptured_groups\tuncaptured_bases\tuncaptured_contigs\tstatus");
			summary.writer.print(b.clear().append(representatives.size()).tab().append(reused).tab().append(additional)
				.tab().append(bases).tab().append(records).append("\tASSEMBLY_SELECTION_NOT_EXACT_SEQUENCE_HOLDOUT_PROOF\n"));
		}
	}

	/** Retained identities can contain UTF-8; LineParser1.parseString is ASCII. */
	private String text(final int field){
		p.setBounds(field);
		final String result=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(result.indexOf('\ufffd')<0, "Invalid UTF-8 would change the assembly contribution identity");
		return result;
	}

	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static final class Sample{
		Sample(final int id_, final String fasta_, final long bases_, final int records_, final byte[] selectionRow_){
			id=id_;
			fasta=fasta_;
			bases=bases_;
			records=records_;
			selectionRow=selectionRow_;
		}
		final int id;
		final String fasta;
		final long bases;
		final int records;
		final byte[] selectionRow;
		int group=-1;
		int capturedSlot=-1;
	}

	private final Sample[] samples;
	private final int[] expected, observed;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ObjectIntMap<String> pathGroups=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> genomeGroups=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class);
	private final ObjectSet<String> samplePaths=new ObjectSet<String>(String.class);
	private final ObjectSet<String> donorGenomes=new ObjectSet<String>(String.class);
	private final IntHashMap sampleIds=new IntHashMap(), groupVoters=new IntHashMap();
	private final IntHashMap uncertainGroups=new IntHashMap(), uncertainDonors=new IntHashMap();
	private int kept, keptExpected, unmappedVoters;
	private boolean selection;
	private static final String[] STATES={"KNOWN_MEMBER_CONTRIBUTOR", "NO_KNOWN_MEMBER_CONTRIBUTOR", "UNRESOLVED_ASSEMBLY_CONTRIBUTION"};
	private static final String ASSIGNMENT_HEADER="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus";
	public static final String SAMPLE_HEADER="input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tactual_attempts\ttarget_loci\tseed_cpu_ns\tnative_index_cpu_ns\toracle_cpu_ns\tranking_cpu_ns\tworker_report_cpu_ns";
	public static final String SELECTION_HEADER="input_id\tfasta\tnucleotide_bases\tfasta_records\tdomain\tphylum\tbase_log2\trecord_log2";
}
