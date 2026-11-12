package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Describes coverage changes from saved native shortlist-sweep tables.
 * Reads metadata only: no FASTA/GFF payload, caller execution or regrading.
 * Input IDs remain distinct; they are not inferred to be independent organisms.
 * Qualified member incidences overlap and must not be summed as distinct genes.
 * One instance owns its readers, retained joins and writers on one thread.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaSweep{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final SummarizeNcrnaSweep report=new SummarizeNcrnaSweep(args);
		try{
			report.process();
			timer.stop();
			report.outstream.println("Summarized "+report.inputs.size()+" original inputs; time: "+timer);
		}finally{
			Shared.closeStream(report.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adapts template.A_SampleByteFile's native CLI and file ownership. */
	public SummarizeNcrnaSweep(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		String requested=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("changes")){
				changes=b;
			}else if(a.equals("summary")){
				summary=b;
			}else if(a.equals("reference")){
				reference=b;
			}else if(a.equals("loci")){
				loci=b;
			}else if(a.equals("members")){
				members=b;
			}else if(a.equals("family")){
				family=b;
			}else if(a.equals("policies")){
				requested=b;
			}else if(a.equals("baseline")){
				baseline=b;
			}else if(a.equals("highlight")){
				highlight=b;
			}else if(a.equals("taxonomy")){
				taxonomy=b;
			}else if(a.equals("compare")){
				require(b!=null && !b.isEmpty(), "compare requires two comma-separated policy names");
				comparison=b;
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		in1=parser.in1;
		prefix=parser.out1;
		overwrite=parser.overwrite;
		require(!parser.append && parser.in2==null && parser.out2==null, "Append and paired files do not apply to saved sweep metadata");
		require(in1!=null && changes!=null && summary!=null && reference!=null && loci!=null && members!=null,
			"Supply in, changes, summary, reference, loci and members metadata files");
		require(prefix!=null && family!=null && !family.isEmpty() && baseline!=null && !baseline.isEmpty()
			&& requested!=null && !requested.isEmpty(), "Supply out=prefix, family and policies; baseline must be named");
		final String[] names=requested.split(",", -1); // One-time CLI list, never a row parser.
		policies=new String[names.length+1];
		policies[0]=baseline;
		System.arraycopy(names, 0, policies, 1, names.length);
		policyTokens=new byte[policies.length][];
		final ObjectIntMap<String> unique=new ObjectIntMap<String>(String.class);
		for(int i=0; i<policies.length; i++){
			require(!policies[i].isEmpty() && !unique.contains(policies[i]), "Policies must be nonempty, distinct and exclude baseline: "+policies[i]);
			unique.put(policies[i], i);
			policyTokens[i]=field(new ByteBuilder(), policies[i]).toBytes();
		}
		if(comparison!=null){
			final String[] pair=comparison.split(",", -1);
			require(taxonomy!=null && pair.length==2 && unique.contains(pair[0]) && unique.contains(pair[1]),
				"compare requires taxonomy and two named policies in baseline/policies");
			leftPolicy=unique.get(pair[0]);
			rightPolicy=unique.get(pair[1]);
			require(leftPolicy!=rightPolicy, "Compared policies must be distinct");
		}
		final String[] paths=new String[6+SUFFIXES.length+(taxonomy==null ? 0 : 2)+(comparison==null ? 0 : 2)];
		paths[0]=in1;
		paths[1]=changes;
		paths[2]=summary;
		paths[3]=reference;
		paths[4]=loci;
		paths[5]=members;
		for(int i=0; i<SUFFIXES.length; i++){paths[6+i]=prefix+SUFFIXES[i];}
		if(taxonomy!=null){
			final int optional=6+SUFFIXES.length;
			paths[optional]=taxonomy;
			paths[optional+1]=prefix+".strata.tsv";
			require(Tools.testInputFiles(false, true, taxonomy), "Cannot read saved taxonomy categories");
			require(Tools.testOutputFiles(overwrite, false, false, paths[optional+1]), "Cannot write category report");
		}
		if(comparison!=null){
			paths[paths.length-2]=prefix+".comparison.loci.tsv";
			paths[paths.length-1]=prefix+".comparison.summary.tsv";
			require(Tools.testOutputFiles(overwrite, false, false, paths[paths.length-2], paths[paths.length-1]),
				"Cannot write comparison reports");
		}
		require(Tools.testInputFiles(false, true, in1, changes, summary, reference, loci, members), "Cannot read metadata inputs");
		require(Tools.testForDuplicateFiles(true, paths), "Input and output paths must be distinct");
		require(Tools.testOutputFiles(overwrite, false, false, paths[6], paths[7], paths[8], paths[9]), "Cannot write report outputs");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/** All joins and conservation checks finish before any output opens. */
	private void process(){
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){
			final Genome genome=new Genome(binding, policies.length);
			inputs.add(genome);
			byId.put(binding.id, genome);
		}
		require(!inputs.isEmpty(), "Reference manifest must contain the sweep population");
		readMembers();
		readCounts();
		readChanges();
		readSupport();
		validate();
		if(taxonomy!=null){readTaxonomy();}
		writeResults();
		if(taxonomy!=null){writeStrata();}
		if(comparison!=null){writeComparison();}
	}

	/** Uses saved full-corpus categories, never recomputes rarity from this panel.
	 * Off-panel rows are ignored. Every on-panel ID must join its exact FASTA once.
	 */
	private void readTaxonomy(){
		assert(taxonomy!=null) : "Taxonomy is optional; the default metadata contract opens no category file";
		try(Table table=new Table(taxonomy)){
			final int id=table.column("input_id"), fasta=table.column("fasta"), category=table.column("category");
			while(table.next()){
				final Genome genome=input(table, id);
				if(genome==null){continue;}
				require(genome.category<0 && genome.binding.fasta.equals(table.text(fasta)),
					"Taxonomy must bind each on-panel input and exact FASTA once: "+genome.binding.id);
				for(int c=0; c<CATEGORIES.length; c++){
					if(table.parser.termEquals(CATEGORIES[c], category)){
						genome.category=c;
						break;
					}
				}
				require(genome.category>=0, "Unknown saved taxonomy category for "+genome.binding.id);
			}
		}
		for(final Genome genome:inputs){require(genome.category>=0, "Missing saved taxonomy input: "+genome.binding.id);}
	}

	/** Counts first and last detections per original input, not changed-copy incidences. */
	private void writeStrata(){
		assert(taxonomy!=null) : "Category output requires a complete validated taxonomy join";
		try(Output out=new Output(FileFormat.testOutput(prefix+".strata.tsv", FileFormat.TXT, null, true, overwrite, false, false))){
			out.writer.println("# Categories are supplied by the full-corpus taxonomy table. UNKNOWN/UNCLASSIFIED remain in aggregate breadth; they are not rare or common classifications.");
			out.writer.println("policy\tcategory\tfirst_gain_inputs\tlast_loss_inputs");
			for(int p=0; p<policies.length; p++){
				final long[][] counts=new long[CATEGORIES.length][2];
				for(final Genome genome:inputs){
					final long base=genome.counts[0][TP], now=genome.counts[p][TP];
					if(base==0 && now>0){counts[genome.category][0]++;}
					if(base>0 && now==0){counts[genome.category][1]++;}
				}
				for(int c=0; c<CATEGORIES.length; c++){
					field(buffer.clear(), policies[p]).tab().append(CATEGORIES[c]).tab().append(counts[c][0]).tab().append(counts[c][1]).nl();
					out.writer.print(buffer);
				}
			}
		}
	}

	/** Joins exact scored-reference keys, never emitted-call or biological-gene sets.
	 * Parents unchanged under every requested policy need not be retained: they
	 * belong to both detection sets or neither and cannot affect set inclusion.
	 */
	private void writeComparison(){
		assert(leftPolicy>=0 && rightPolicy>=0 && taxonomy!=null)
			: "Comparison uses resolved policy indices and the validated original-input category join";
		final long[][] counts=new long[CATEGORIES.length+1][8];
		try(Output lociOut=comparisonOutput("loci"); Output summaryOut=comparisonOutput("summary")){
			lociOut.writer.println("# Exact keys are original input_id plus parent_id. Rows include any reference changed by either compared policy; unchanged references are omitted. Coordinates are 1-based inclusive.");
			lociOut.writer.println("input_id\tparent_id\tcategory\tfamily\tseqid\tstart\tstop\tstrand\tbaseline_present\tleft_policy\tleft_change\tleft_present\tright_policy\tright_change\tright_present");
			summaryOut.writer.println("# Detection sets contain scored-reference keys, not emitted GFF calls. Categories are the supplied full-corpus labels; ALL sums original input records without biological deduplication.");
			summaryOut.writer.println("left_policy\tright_policy\tcategory\tleft_TP\tright_TP\tleft_lost\tleft_gained\tright_lost\tright_gained\tleft_only\tright_only\tleft_subset_right\tright_subset_left");
			for(final Genome genome:inputs){
				final long[] row=counts[genome.category];
				for(int side=0; side<2; side++){
					final long[] nativeCounts=genome.counts[side==0 ? leftPolicy : rightPolicy];
					row[side]=Math.addExact(row[side], nativeCounts[TP]);
					row[2+2*side]=Math.addExact(row[2+2*side], nativeCounts[LOST]);
					row[3+2*side]=Math.addExact(row[3+2*side], nativeCounts[GAINED]);
				}
				for(final Parent parent:genome.parents){
					if(parent.directions[leftPolicy]==0 && parent.directions[rightPolicy]==0){continue;}
					final boolean left=present(parent, leftPolicy), right=present(parent, rightPolicy);
					if(left && !right){row[6]++;}
					if(right && !left){row[7]++;}
					writeComparisonLocus(lociOut, genome, parent, left, right);
				}
			}
			final long[] all=counts[CATEGORIES.length];
			for(int c=0; c<CATEGORIES.length; c++){
				for(int i=0; i<all.length; i++){all[i]=Math.addExact(all[i], counts[c][i]);}
			}
			for(int c=0; c<counts.length; c++){
				final long[] row=counts[c];
				require(row[1]-row[0]==row[7]-row[6],
					"Exact reference-set differences must conserve native TP totals for comparison category "+c);
				field(buffer.clear(), policies[leftPolicy]).tab();
				field(buffer, policies[rightPolicy]).tab().append(c<CATEGORIES.length ? CATEGORIES[c] : "ALL");
				for(final long count:row){buffer.tab().append(count);}
				summaryOut.writer.print(buffer.tab().append(row[6]==0).tab().append(row[7]==0).nl());
			}
		}
	}

	/** Baseline membership is fixed by the single validated loss/gain direction. */
	private static boolean present(final Parent parent, final int policy){
		assert(parent.direction==LOST || parent.direction==GAINED)
			: "readChanges retains only changed parents and rejects opposing directions relative to one baseline";
		return parent.directions[policy]==0 ? parent.direction==LOST : parent.directions[policy]==GAINED;
	}

	private void writeComparisonLocus(final Output out, final Genome genome, final Parent parent,
			final boolean left, final boolean right){
		field(buffer.clear(), genome.binding.id).tab();
		field(buffer, parent.id).tab().append(CATEGORIES[genome.category]).tab();
		field(buffer, family).tab();
		field(buffer, parent.seqid).tab().append(parent.start).tab().append(parent.stop).tab().append(parent.strand)
			.tab().append(parent.direction==LOST).tab();
		field(buffer, policies[leftPolicy]).tab().append(changeName(parent.directions[leftPolicy])).tab().append(left).tab();
		field(buffer, policies[rightPolicy]).tab().append(changeName(parent.directions[rightPolicy])).tab().append(right).nl();
		out.writer.print(buffer);
	}

	private static String changeName(final int direction){
		return direction==LOST ? "LOST" : direction==GAINED ? "GAINED" : "UNCHANGED";
	}

	private Output comparisonOutput(final String kind){
		return new Output(FileFormat.testOutput(prefix+".comparison."+kind+".tsv", FileFormat.TXT, null, true, overwrite, false, false));
	}

	/** The member map is literal UTF-8, matching NcrnaSubtypeEvidence's contract. */
	private void readMembers(){
		try(Table table=new Table(members)){
			final int f=table.column("family"), m=table.column("member"), r=table.column("reference_model");
			while(table.next()){
				if(!table.literal(f).equals(family)){continue;}
				final String name=table.literal(m);
				require(!name.isEmpty() && !name.equals(".") && name.indexOf(',')<0 && !memberIndex.contains(name),
					"Member names must be distinct and usable in native comma lists: "+name);
				memberIndex.put(name, memberList.size());
				memberList.add(new Member(name, table.literal(r)));
			}
		}
		require(!memberList.isEmpty(), "No declared reference members for family "+family);
		if(highlight!=null){
			for(final String name:highlight.split(",", -1)){
				require(memberIndex.contains(name), "Highlighted member is not declared: "+name);
				memberList.get(memberIndex.get(name)).highlight=true;
			}
		}
	}

	/** Retains only three native integer metrics per requested input/policy. */
	private void readCounts(){
		try(Table table=new Table(in1)){
			final int id=table.column("input_id"), p=table.column("policy");
			final int tp=table.column("TP"), lost=table.column("lost_baseline_TP"), gained=table.column("gained_TP");
			while(table.next()){
				final int policy=policy(table, p);
				if(policy<0){continue;}
				final Genome genome=input(table, id);
				require(genome!=null, "Sweep count has an input outside the reference manifest");
				if(genome.seen[policy]){
					throw new IllegalArgumentException("Repeated input/policy count: "+genome.binding.id+" / "+policies[policy]);
				}
				genome.seen[policy]=true;
				genome.counts[policy][TP]=table.count(tp);
				genome.counts[policy][LOST]=table.count(lost);
				genome.counts[policy][GAINED]=table.count(gained);
			}
		}
	}

	/** A parent is retained once per original input, with a direction per policy. */
	private void readChanges(){
		try(Table table=new Table(changes)){
			final int id=table.column("input_id"), p=table.column("policy"), change=table.column("change");
			final int ref=table.column("reference_id"), seq=table.column("seqid");
			final int start=table.column("start"), stop=table.column("stop"), strand=table.column("strand");
			while(table.next()){
				final int policy=policy(table, p);
				if(policy<0){continue;}
				require(policy>0, "Baseline cannot contain paired changes against itself");
				final Genome genome=input(table, id);
				require(genome!=null, "Changed parent has an input outside the reference manifest");
				final byte direction=(byte)(table.parser.termEquals("LOST", change) ? LOST
					: table.parser.termEquals("GAINED", change) ? GAINED : 0);
				require(direction!=0, "Change must be LOST or GAINED");
				final String parentId=table.text(ref);
				Parent parent=genome.byParent.get(parentId);
				if(parent==null){
					parent=new Parent(parentId, table.text(seq), table.parser.parseInt(start),
						table.parser.parseInt(stop), table.text(strand), policies.length);
					genome.byParent.put(parentId, parent);
					genome.parents.add(parent);
				}else{
					geometry(parent, table, seq, start, stop, strand);
				}
				require(parent.directions[policy]==0, "Repeated changed parent/policy: "+genome.binding.id+" / "+parentId+" / "+policies[policy]);
				require(parent.direction==0 || parent.direction==direction,
					"A parent cannot be both LOST and GAINED relative to one baseline: "+genome.binding.id+" / "+parentId);
				parent.direction=direction;
				parent.directions[policy]=direction;
				genome.counts[policy][direction+2]++;
			}
		}
	}

	/** Streams the full inventory but retains support only for changed parents. */
	private void readSupport(){
		try(Table table=new Table(loci)){
			final int id=table.column("input_id"), ref=table.column("parent_id"), f=table.column("family");
			final int seq=table.column("seqid"), start=table.column("start"), stop=table.column("stop"), strand=table.column("strand");
			final int status=table.column("status"), scored=table.column("scored_members"), marginal=table.column("marginal_members");
			final byte[] familyToken=field(new ByteBuilder(), family).toBytes();
			while(table.next()){
				if(!table.matches(familyToken, f)){continue;}
				final Genome genome=input(table, id);
				if(genome==null || genome.parents.isEmpty()){continue;}
				Parent parent=null;
				for(final Parent candidate:genome.parents){
					if(table.matches(candidate.token, ref)){
						parent=candidate;
						break;
					}
				}
				if(parent==null){continue;}
				require(parent.scored==null, "Duplicate inventory binding for changed parent: "+parent.id);
				require(table.parser.termEquals("scored", status), "Changed TP must bind to a scored inventory parent: "+parent.id);
				geometry(parent, table, seq, start, stop, strand);
				parent.scored=memberSupport(table, scored);
				parent.marginal=memberSupport(table, marginal);
				boolean qualified=false;
				for(int m=0; m<memberList.size(); m++){
					qualified|=parent.scored[m];
					require(!parent.scored[m] || !parent.marginal[m], "A member cannot be both qualified and marginal for the same parent: "+parent.id);
				}
				require(qualified, "Scored parent requires qualified member support, as enforced by NcrnaSubtypeEvidence.support: "+parent.id);
			}
		}
	}

	/** Decodes member tokens separately; marginal evidence never becomes scored. */
	private boolean[] memberSupport(final Table table, final int column){
		final boolean[] support=new boolean[memberList.size()];
		if(table.parser.termEquals('.', column)){return support;}
		table.parser.setBounds(column);
		final byte[] line=table.parser.line();
		final int end=table.parser.b();
		for(int begin=table.parser.a(); begin<end;){
			int stop=begin;
			while(stop<end && line[stop]!=','){stop++;}
			final String name=table.decoder.decode(line, begin, stop);
			require(memberIndex.contains(name), "Inventory support names an undeclared member: "+name);
			final int index=memberIndex.get(name);
			require(!support[index], "Repeated member within one support list: "+name);
			support[index]=true;
			require(stop==end || stop+1<end, "Trailing empty inventory member");
			begin=stop+1;
		}
		return support;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Accounting          ----------------*/
	/*--------------------------------------------------------------*/

	/** Cross-checks paired changes against both per-input and whole-policy totals. */
	private void validate(){
		totals=new long[policies.length][3];
		for(final Genome genome:inputs){
			for(int p=0; p<policies.length; p++){
				if(!genome.seen[p]){
					throw new IllegalArgumentException("Missing input/policy count: "+genome.binding.id+" / "+policies[p]);
				}
				final long[] row=genome.counts[p];
				if(row[LOST]!=row[3] || row[GAINED]!=row[4]){
					throw new IllegalArgumentException("Paired changes disagree with native counts for "+genome.binding.id+" / "+policies[p]);
				}
				if(row[LOST]>genome.counts[0][TP] || row[GAINED]>row[TP]
					|| Math.addExact(genome.counts[0][TP]-row[LOST], row[GAINED])!=row[TP]){
					throw new IllegalArgumentException("Baseline TP minus lost plus gained must equal policy TP for "+genome.binding.id+" / "+policies[p]);
				}
				for(int c=0; c<3; c++){totals[p][c]=Math.addExact(totals[p][c], row[c]);}
			}
			for(final Parent parent:genome.parents){
				require(parent.scored!=null, "No scored inventory binding for "+genome.binding.id+" / "+parent.id);
			}
		}
		final boolean[] seen=new boolean[policies.length];
		try(Table table=new Table(summary)){
			final int p=table.column("policy"), tp=table.column("TP");
			final int lost=table.column("lost_baseline_TP"), gained=table.column("gained_TP");
			while(table.next()){
				final int policy=policy(table, p);
				if(policy<0){continue;}
				require(!seen[policy], "Repeated policy summary: "+policies[policy]);
				seen[policy]=true;
				require(totals[policy][TP]==table.count(tp) && totals[policy][LOST]==table.count(lost)
					&& totals[policy][GAINED]==table.count(gained), "Input and paired counts disagree with policy summary: "+policies[policy]);
			}
		}
		for(int p=0; p<seen.length; p++){require(seen[p], "Missing policy summary: "+policies[p]);}
	}

	/** Maps one incidence to a detection-coverage class, not sequence redundancy. */
	private static int coverage(final long baselineTP, final long policyTP, final int direction){
		assert(direction==LOST || direction==GAINED) : "Only paired loss/gain incidences receive coverage labels";
		return direction==LOST ? (policyTP==0 ? 1 : 0) : (baselineTP==0 ? 2 : 3);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits all reports in manifest order and requested-policy order. */
	private void writeResults(){
		try(Output locusOut=output(0); Output inputOut=output(1); Output summaryOut=output(2); Output memberOut=output(3)){
			locusOut.writer.println("policy\tchange\tinput_id\tparent_id\tfamily\tseqid\tstart\tstop\tstrand\tbaseline_TP\tpolicy_TP\tcoverage\tfasta\tscored_members\tmarginal_members");
			inputOut.writer.println("policy\tinput_id\tbaseline_TP\tpolicy_TP\tlost_loci\tgained_loci\tlost_coverage\trecovered_coverage\tfasta");
			summaryOut.writer.println("policy\tinputs\tbaseline_TP\tpolicy_TP\tlost_loci\tgained_loci\taffected_inputs\tretained_loss_loci\tretained_loss_inputs\tlast_lost_loci\tlast_lost_inputs\tfirst_recovered_loci\tfirst_recovered_inputs\tadditional_gain_loci\tadditional_gain_inputs");
			memberOut.writer.println("policy\tchange\tmember\treference_model\thighlight\tscored_loci\tscored_inputs\tmarginal_loci\tmarginal_inputs\tqualification");
			for(int p=1; p<policies.length; p++){
				final long[] classLoci=new long[4], classInputs=new long[4];
				final long[][][] memberCounts=new long[2][memberList.size()][4];
				final boolean[][][] memberSeen=new boolean[2][memberList.size()][2];
				long affected=0;
				for(final Genome genome:inputs){
					final long[] row=genome.counts[p];
					if(row[LOST]+row[GAINED]==0){continue;}
					affected++;
					final long base=genome.counts[0][TP];
					writeInput(inputOut, genome, p);
					for(int direction=LOST; direction<=GAINED; direction++){
						if(row[direction]==0){continue;}
						final int c=coverage(base, row[TP], direction);
						classLoci[c]+=row[direction];
						classInputs[c]++;
					}
					for(final boolean[][] direction:memberSeen){
						for(final boolean[] member:direction){Arrays.fill(member, false);}
					}
					for(final Parent parent:genome.parents){
						final int direction=parent.directions[p];
						if(direction==0){continue;}
						writeLocus(locusOut, genome, parent, p);
						for(int m=0; m<memberList.size(); m++){
							final long[] counts=memberCounts[direction-1][m];
							final boolean[] seen=memberSeen[direction-1][m];
							if(parent.scored[m]){
								counts[0]++;
								seen[0]=true;
							}
							if(parent.marginal[m]){
								counts[2]++;
								seen[1]=true;
							}
						}
					}
					for(int d=0; d<2; d++){
						for(int m=0; m<memberList.size(); m++){
							if(memberSeen[d][m][0]){memberCounts[d][m][1]++;}
							if(memberSeen[d][m][1]){memberCounts[d][m][3]++;}
						}
					}
				}
				assert(classLoci[0]+classLoci[1]==totals[p][LOST] && classLoci[2]+classLoci[3]==totals[p][GAINED])
					: "Coverage classes must partition the validated paired incidences";
				field(buffer.clear(), policies[p]).tab().append(inputs.size()).tab().append(totals[0][TP]);
				for(int c=0; c<3; c++){buffer.tab().append(totals[p][c]);}
				buffer.tab().append(affected);
				for(int c=0; c<4; c++){buffer.tab().append(classLoci[c]).tab().append(classInputs[c]);}
				summaryOut.writer.print(buffer.nl());
				writeMembers(memberOut, p, memberCounts);
			}
		}
	}

	private void writeInput(final Output out, final Genome genome, final int policy){
		final long base=genome.counts[0][TP];
		final long[] row=genome.counts[policy];
		field(buffer.clear(), policies[policy]).tab();
		field(buffer, genome.binding.id).tab().append(base).tab().append(row[TP]).tab().append(row[LOST]).tab().append(row[GAINED])
			.tab().append(base>0 && row[TP]==0).tab().append(base==0 && row[TP]>0).tab();
		field(buffer, genome.binding.fasta).nl();
		out.writer.print(buffer);
	}

	private void writeLocus(final Output out, final Genome genome, final Parent parent, final int policy){
		final int direction=parent.directions[policy];
		final long base=genome.counts[0][TP], tp=genome.counts[policy][TP];
		field(buffer.clear(), policies[policy]).tab().append(direction==LOST ? "LOST" : "GAINED").tab();
		field(buffer, genome.binding.id).tab();
		field(buffer, parent.id).tab();
		field(buffer, family).tab();
		field(buffer, parent.seqid).tab().append(parent.start).tab().append(parent.stop).tab().append(parent.strand)
			.tab().append(base).tab().append(tp).tab().append(CLASSES[coverage(base, tp, direction)]).tab();
		field(buffer, genome.binding.fasta).tab();
		appendMembers(parent.scored);
		buffer.tab();
		appendMembers(parent.marginal);
		out.writer.print(buffer.nl());
	}

	private void writeMembers(final Output out, final int policy, final long[][][] counts){
		for(int d=0; d<2; d++){
			for(int m=0; m<memberList.size(); m++){
				final Member member=memberList.get(m);
				field(buffer.clear(), policies[policy]).tab().append(d==0 ? "LOST" : "GAINED").tab();
				field(buffer, member.name).tab();
				field(buffer, member.model).tab().append(member.highlight);
				for(final long count:counts[d][m]){buffer.tab().append(count);}
				out.writer.print(buffer.append("\tNONADDITIVE_REFERENCE_SUPPORT\n"));
			}
		}
	}

	private void appendMembers(final boolean[] support){
		boolean first=true;
		for(int m=0; m<support.length; m++){
			if(!support[m]){continue;}
			if(!first){buffer.append(',');}
			field(buffer, memberList.get(m).name);
			first=false;
		}
		if(first){buffer.append('.');}
	}

	private Output output(final int index){
		return new Output(FileFormat.testOutput(prefix+SUFFIXES[index], FileFormat.TXT, null, true, overwrite, false, false));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Shared Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private int policy(final Table table, final int column){
		for(int i=0; i<policies.length; i++){
			if(table.matches(policyTokens[i], column)){return i;}
		}
		return -1;
	}

	/** Sweep/inventory rows group inputs; cache the last token, including misses. */
	private Genome input(final Table table, final int column){
		if(table.lastInput==null || !table.matches(table.lastInput, column)){
			table.lastInput=table.parser.parseByteArray(column);
			table.genome=byId.get(table.text(column));
		}
		return table.genome;
	}

	private static void geometry(final Parent parent, final Table table, final int seq, final int start, final int stop, final int strand){
		require(parent.seqid.equals(table.text(seq)) && parent.start==table.parser.parseInt(start)
			&& parent.stop==table.parser.parseInt(stop) && parent.strand.equals(table.text(strand)),
			"Changed parent geometry disagrees across metadata: "+parent.id);
	}
	private static ByteBuilder field(final ByteBuilder out, final String text){return BuildNcrnaReference.field(out, text);}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Row Owners          ----------------*/
	/*--------------------------------------------------------------*/

	/** Header-indexed native byte reader; only retained text fields allocate. */
	static final class Table implements AutoCloseable{
		Table(final String path_){
			path=path_;
			input=new Input(path);
			try{
				byte[] line;
				do{line=input.reader.nextLine();}while(line!=null && (line.length==0 || line[0]=='#'));
				require(line!=null, "Missing table header: "+path);
				parser.set(line);
				width=parser.terms();
				for(int i=0; i<width; i++){
					final String name=literal(i);
					require(!columns.contains(name), "Duplicate header column: "+name+" in "+path);
					columns.put(name, i);
				}
			}catch(RuntimeException | Error failure){
				try{input.close();}catch(RuntimeException | Error close){failure.addSuppressed(close);}
				throw failure;
			}
		}
		boolean next(){
			byte[] line;
			do{line=input.reader.nextLine();}while(line!=null && (line.length==0 || line[0]=='#'));
			if(line==null){return false;}
			parser.set(line);
			if(parser.terms()!=width){
				throw new IllegalArgumentException("Wrong table width in "+path+" at line "+input.reader.lineNum());
			}
			return true;
		}
		int column(final String name){
			require(columns.contains(name), "Missing column "+name+" in "+path);
			return columns.get(name);
		}
		String literal(final int column){
			parser.setBounds(column);
			return new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		}
		String text(final int column){
			parser.setBounds(column);
			return decoder.decode(parser.line(), parser.a(), parser.b());
		}
		/** Byte comparison preserves UTF-8 identity without a temporary field String. */
		boolean matches(final byte[] token, final int column){
			parser.setBounds(column);
			if(parser.b()-parser.a()!=token.length){return false;}
			final byte[] line=parser.line();
			for(int i=0; i<token.length; i++){
				if(token[i]!=line[parser.a()+i]){return false;}
			}
			return true;
		}
		long count(final int column){
			final long count=parser.parseLong(column);
			if(count<0){throw new IllegalArgumentException("Negative count in "+path);}
			return count;
		}
		@Override
		public void close(){input.close();}
		final String path;
		final Input input;
		final LineParser1 parser=new LineParser1('\t');
		final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
		final ObjectIntMap<String> columns=new ObjectIntMap<String>(String.class);
		int width;
		byte[] lastInput;
		Genome genome;
	}

	private static final class Genome{
		Genome(final Binding binding_, final int policies){
			binding=binding_;
			counts=new long[policies][5];
			seen=new boolean[policies];
		}
		final Binding binding;
		/** TP, native lost/gained, then observed paired lost/gained. */
		final long[][] counts;
		final boolean[] seen;
		final ObjectMap<String, Parent> byParent=new ObjectMap<String, Parent>(String.class, Parent.class);
		final ArrayList<Parent> parents=new ArrayList<Parent>();
		/** Index in the supplied full-corpus category table; -1 until joined. */
		int category=-1;
	}
	private static final class Parent{
		Parent(final String id_, final String seqid_, final int start_, final int stop_, final String strand_, final int policies){
			id=id_;
			token=field(new ByteBuilder(), id).toBytes();
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			require(!id.isEmpty() && !seqid.isEmpty() && start>=1 && stop>=start && (strand.equals("+") || strand.equals("-")),
				"Parents require IDs, 1-based inclusive coordinates and a strand: "+id);
			directions=new byte[policies];
		}
		final String id, seqid, strand;
		final byte[] token;
		final int start, stop;
		final byte[] directions;
		/** Baseline membership fixes the direction across all compared policies. */
		byte direction;
		boolean[] scored, marginal;
	}
	private static final class Member{
		Member(final String name_, final String model_){
			name=name_;
			model=model_;
		}
		final String name, model;
		boolean highlight;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1, changes, summary, reference, loci, members, family, prefix, highlight, taxonomy;
	private String baseline="baseline";
	/** Optional exact detection-set comparison; both indices resolve before I/O. */
	private String comparison;
	private int leftPolicy=-1, rightPolicy=-1;
	private final String[] policies;
	private final byte[][] policyTokens;
	private final boolean overwrite;
	private final PrintStream outstream;
	private final ObjectMap<String, Genome> byId=new ObjectMap<String, Genome>(String.class, Genome.class);
	private final ArrayList<Genome> inputs=new ArrayList<Genome>();
	private final ArrayList<Member> memberList=new ArrayList<Member>();
	private final ObjectIntMap<String> memberIndex=new ObjectIntMap<String>(String.class);
	private final ByteBuilder buffer=new ByteBuilder();
	private long[][] totals;
	private static final int TP=0, LOST=1, GAINED=2;
	private static final String[] SUFFIXES={".loci.tsv", ".inputs.tsv", ".summary.tsv", ".members.tsv"};
	private static final String[] CATEGORIES={"CORPUS_LOWER_QUARTILE", "OTHER_CLASSIFIED", "UNKNOWN", "UNCLASSIFIED"};
	private static final String[] CLASSES={"OTHER_TARGET_DETECTIONS_REMAIN", "LAST_DETECTED_COPY_LOST",
		"FIRST_DETECTION_RECOVERED", "ADDITIONAL_DETECTION_IN_COVERED_INPUT"};
}
