package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectIntMap;
import prok.NcrnaConstructionMembers;
import prok.BuildNcrnaDdlProfiles;
import prok.NcrnaDdlHoldout;
import prok.NcrnaDdlMembers;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaExactSequenceSet;
import prok.NcrnaLearningExclusions;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;

/** Current-roster alias, literal-IUPAC, removed-cluster and empty-profile fixture.
 * No real genome data or source/TaxID inference enters this tiny population.
 * @author Brian, Raiden, Ganyu
 */
public final class NcrnaDdlHoldoutTest{

	public static void run(final String root){
		check(new File(root).mkdir(), "The holdout fixture directory must be fresh");
		final String manifest=root+"/manifest.tsv", panel=root+"/panel.tsv", sources=root+"/sources.tsv", aliases=root+"/aliases.tsv";
		final String fasta=root+"/members.fa", assignments=root+"/assignments.tsv";
		final String[] ids={"RF00174_input_A", "RF00174_input_B", "RF00174_input_C", "RF00174_input_D", "RF00174_input_E"};
		final String[] dna={"ACGRNT", "GATCGATC", "acgrnt", "TTTACG", "CTAGGA"};
		final String[] names={"cobalamin_RF00174_zero", "cobalamin_RF00174_retained"};
		try(Output out=new Output(manifest)){
			out.writer.println(NcrnaDdlHoldout.MANIFEST_HEADER);
			for(char letter='A'; letter<='D'; letter++){
				out.writer.println("assembly_"+letter+"\t"+letter+"\tFIXTURE\t"+root+"/"+letter+".fa\t"+root+"/"+letter+
					".gff\t.\tACC_"+letter+"\t"+root+"/source_"+letter+".fa\t.\tFASTA\tNCRNA\tCORE");
			}
		}
		try(Output out=new Output(aliases)){
			out.writer.println(NcrnaConstructionMembers.GROUP_HEADER);
			for(char letter='A'; letter<='D'; letter++){
				final int group=letter<='B' ? 10 : letter=='C' ? 20 : 30;
				out.writer.println("FIXTURE\t"+letter+"\t"+root+"/source_"+letter+".fa\t"+root+"/"+letter+
					".fa\tACC_"+letter+"\t"+group+"\t0\tfalse\tfalse\tRESOLVED\t"+group+"\t0\tfalse\tRESOLVED");
			}
		}
		try(Output out=new Output(panel)){
			out.writer.println(NcrnaDdlHoldout.PANEL_HEADER);
			out.writer.println("0\tassembly_A\tA\t"+root+"/A.fa\t"+root+"/A.gff");
		}
		try(Output out=new Output(sources)){
			out.writer.println(NcrnaDdlHoldout.SOURCE_HEADER);
			for(int i=0; i<ids.length; i++){
				final char genome=(char)('A'+Math.min(i, 3));
				out.writer.println(ids[i]+"\tRF00174\tassembly_"+genome+"\t"+genome+
					"\tcontig\t1\t"+dna[i].length()+"\t+\t"+dna[i].length()+"\t.\t0");
			}
		}
		try(Output out=new Output(assignments)){
			out.writer.println("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus");
			for(int i=0; i<ids.length; i++){
				final String model=i<2 ? names[0] : i<4 ? names[1] : "cobalamin_RF00174_removed";
				out.writer.println(ids[i]+" details=kept\t0\tRF00174\tRF00174_"+i+"\t"+model+"\tkept");
			}
		}
		writeFasta(fasta, ids, dna);
		literal();
		final ObjectIntMap<String> byAssembly=new ObjectIntMap<String>(String.class), byPath=new ObjectIntMap<String>(String.class);
		final ArrayList<String> assemblies=new ArrayList<String>();
		final int[] rows={0};
		NcrnaLearningExclusions.readManifest(manifest, byAssembly, byPath, assemblies, (id, linked, source, row)->rows[0]++);
		check(rows[0]==4 && assemblies.size()==4 && byAssembly.get("assembly_D")==3
			&& byPath.get(root+"/A.fa")==0 && byPath.get(root+"/source_A.fa")==0,
			"The shared manifest parser must retain the existing exact assembly/path indexes and row order");
		final NcrnaDdlHoldout holdout=new NcrnaDdlHoldout(manifest, panel, sources, aliases);
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, iupac=Read.IUPAC_TO_N;
		try{
			Shared.TRIM_READ_DESCRIPTION=false;
			Read.IUPAC_TO_N=false;
			final byte[][] references={bytes("ACGRNT"), bytes("TTTACG")};
			final int[] capturedEdges={0};
			final NcrnaDdlMembers members=NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, null, names, references, 5, holdout,
				(id, model, sequence)->{
					check(id.startsWith(ids[3]+" ") && model==1 && java.util.Arrays.equals(sequence, bytes(dna[3])),
						"Training incidence observes only the exact retained original/model edge, never a held-out or removed assignment");
					capturedEdges[0]++;
				});
			check(capturedEdges[0]==1 && holdout.allowedOriginal(ids[3], bytes(dna[3]))
				&& !holdout.allowedOriginal(ids[0], bytes(dna[0])) && !holdout.allowedOriginal(ids[2], bytes(dna[2]))
				&& holdout.allowedTrainingAssembly("assembly_D") && !holdout.allowedTrainingAssembly("assembly_A")
				&& !holdout.allowedTrainingAssembly("assembly_B") && !holdout.allowedTrainingAssembly("assembly_C"),
				"Read-only GC eligibility retains panel, alias and exact-sequence exclusions without changing voter counters");
			check(members.inputMembers()==5 && members.keptMembers()==1 && members.excludedMembers()==3
				&& members.removedMembers()==1 && members.members(0).length==0 && members.members(1).length==1
				&& members.assignedMembers(0)==2 && members.assignedMembers(1)==2,
				"Panel/alias and exact-DNA exclusions must conserve originals, retain an empty model and never reassign a removed cluster");
			refitBuilder(root, names, references, dna, manifest, panel, sources, aliases, fasta, assignments);
			final String[] unbound={names[0], names[1], "cobalamin_RF00174_unbound"};
			reject(()->NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, unbound,
				new byte[][]{references[0], references[1], references[0]}, 5,
				new NcrnaDdlHoldout(manifest, panel, sources, aliases)),
				"A missing current-model assignment must fail rather than masquerade as a holdout-empty model");
			final NcrnaDdlHoldout repeated=new NcrnaDdlHoldout(manifest, panel, sources, aliases);
			repeated.observe(ids[0]+" header=one", bytes(dna[0]));
			reject(()->repeated.observe(ids[0]+" header=two", bytes(dna[0])),
				"Distinct full headers may not repeat a first-token original source identity");
			final NcrnaDdlHoldout incomplete=new NcrnaDdlHoldout(manifest, panel, sources, aliases);
			incomplete.observe(ids[0], bytes(dna[0]));
			reject(incomplete::finishObservation, "An incomplete main source ledger must fail before voter exclusion");
			holdout.writeSample(root+"/training.tsv", 1, 20261119);
			holdout.writeAudit(root+"/audit.tsv");
			try(Input in=new Input(root+"/training.tsv")){
				check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).equals(NcrnaDdlHoldout.PANEL_HEADER),
					"Training inputs retain the ordinary source-qualified panel schema");
				check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).startsWith("0\tassembly_D\tD\t")
					&& in.reader.nextLine()==null,
					"Sample only after assembly/alias/duplicate exclusions; no matching gene or assembly may leak because of FASTA order");
			}
			final String duplicateRemoved=root+"/duplicate_removed.fa";
			final String[] removedDna=dna.clone();
			removedDna[4]=dna[0];
			writeFasta(duplicateRemoved, ids, removedDna);
			final NcrnaDdlHoldout removedHoldout=new NcrnaDdlHoldout(manifest, panel, sources, aliases);
			final NcrnaDdlMembers retained=NcrnaDdlMembers.loadSubsetRoster(duplicateRemoved, assignments, names, references, 5, removedHoldout);
			check(retained.keptMembers()==1 && retained.members(1).length==1,
				"A removed-cluster duplicate must not erase a distinct retained member from the same genome");
			reject(()->removedHoldout.writeSample(root+"/no_training.tsv", 1, 20261119),
				"Padding sampling conservatively excludes a group even when its only panel-DNA duplicate belongs to a removed cluster");
			minor(root, manifest, panel, aliases);
			protectedProfiles(root, manifest, panel, aliases);
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.IUPAC_TO_N=iupac;
		}
		final String conflicting=root+"/conflicting_aliases.tsv";
		try(Output out=new Output(conflicting)){
			out.writer.println(NcrnaConstructionMembers.GROUP_HEADER);
			out.writer.println("FIXTURE\tX\t"+root+"/source_A.fa\t"+root+"/other_A.fa\t.\t10\t0\tfalse\tfalse\tRESOLVED\t10\t0\tfalse\tRESOLVED");
			out.writer.println("FIXTURE\tY\t"+root+"/other_source.fa\t"+root+"/A.fa\t.\t11\t0\tfalse\tfalse\tRESOLVED\t11\t0\tfalse\tRESOLVED");
		}
		reject(()->new NcrnaDdlHoldout(manifest, panel, sources, conflicting),
			"Known source/canonical groups must not silently discard a contradictory alias edge");
		bonusBuilder(root+"/bonus_builder");
	}

	/** The release switch restores exactly A (panel), B (alias) and C (literal
	 * duplicate), while D stays retained and E's removed cluster stays nonvoting. */
	private static void refitBuilder(final String root, final String[] names, final byte[][] references,
			final String[] dna, final String manifest, final String panel, final String sources,
			final String aliases, final String fasta, final String assignments){
		final NcrnaDdlHoldout all=new NcrnaDdlHoldout(manifest, null, sources, aliases, "RF00174", null, 0, true);
		final NcrnaDdlMembers members=NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, names, references, 5, all);
		check(all.panelGenomes()==0 && members.keptMembers()==4 && members.excludedMembers()==0
			&& members.removedMembers()==1 && members.assignedMembers(0)==2 && members.assignedMembers(1)==2,
			"Release refitting restores the three excluded originals without assigning the removed cluster or changing the model axis");
		check(Arrays.equals(members.members(0)[0], bytes(dna[0])) && Arrays.equals(members.members(0)[1], bytes(dna[1]))
			&& Arrays.equals(members.members(1)[0], bytes(dna[2])) && Arrays.equals(members.members(1)[1], bytes(dna[3])),
			"The literal original sequences and original per-model order survive refitting, including the cross-assembly duplicate");
		final String ref=root+"/refit.models.fa";
		try(Output out=new Output(ref)){
			for(int i=0; i<names.length; i++){
				out.writer.println('>'+names[i]);
				out.writer.println(new String(references[i], StandardCharsets.US_ASCII));
			}
		}
		final boolean amino=Shared.AMINO_IN, uracil=Read.U_TO_T, dots=Read.DOT_DASH_X_TO_N, lower=Read.LOWER_CASE_TO_N;
		try{
			Shared.AMINO_IN=false;
			Read.U_TO_T=false;
			Read.DOT_DASH_X_TO_N=false;
			Read.LOWER_CASE_TO_N=false;
			final ArrayList<String> common=new ArrayList<String>(Arrays.asList("orderonly=t", "in="+fasta,
				"assignments="+assignments, "ref="+ref, "manifest="+manifest, "sources="+sources,
				"aliases="+aliases, "expectedinput=5", "expectedmodels=2", "partitionbins=1024",
				"paddingblacklist=none", "seed=12345", "t=1", "ow=f"));
			final String plain=root+"/default.tsv", explicit=root+"/explicit_false.tsv", refit=root+"/all_data.tsv";
			final ArrayList<String> args=new ArrayList<String>(common);
			args.addAll(Arrays.asList("panel="+panel, "expectedpanel=1", "expectedmembers=1", "out="+plain));
			BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));
			args.set(args.size()-1, "out="+explicit);
			args.add("refit=f");
			BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));
			check(sameBytes(plain, explicit), "Omitted refit and refit=f preserve identical member-only profile bytes");
			args.clear();
			args.addAll(common);
			args.addAll(Arrays.asList("refit=t", "expectedpanel=0", "expectedmembers=4", "out="+refit));
			BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));
			final NcrnaDdlProfileResource before=NcrnaDdlProfileResource.load(plain, names, references);
			final NcrnaDdlProfileResource after=NcrnaDdlProfileResource.load(refit, names, references);
			check(before.members(0)==0 && before.members(1)==1 && after.members(0)==2 && after.members(1)==2
				&& after.seed()==before.seed() && after.buckets()==before.buckets(),
				"The actual CLI removes the positive-panel requirement and restores only the excluded voters in the same profile geometry");
			check(new File(refit+".refit.tsv").isFile() && new File(refit+".dispositions.tsv").isFile()
				&& !new File(refit+".holdout.tsv").exists(), "All-data provenance is labelled refit, never a fabricated holdout");
			args.set(args.size()-1, "out="+root+"/invalid_refit.tsv");
			args.add("panel="+panel);
			reject(()->BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()])),
				"A release refit must reject a supplied evaluation panel instead of silently ignoring it");
			args.remove(args.size()-1);
			args.add("prepareonly=t");
			reject(()->BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()])),
				"All-data refitting cannot pretend to prepare a held-out training sample");
		}finally{
			Shared.AMINO_IN=amino;
			Read.U_TO_T=uracil;
			Read.DOT_DASH_X_TO_N=dots;
			Read.LOWER_CASE_TO_N=lower;
		}
	}

	/** The actual CLI must preserve a held-out zero-voter model and apply its
	 * declared bonus inside the profile; defaults remain byte-identical. */
	private static void bonusBuilder(final String root){
		check(new File(root).mkdir(), "The builder integration fixture needs a fresh directory");
		final String manifest=root+"/manifest.tsv", panel=root+"/panel.tsv";
		final String sources=root+"/sources.tsv", aliases=root+"/aliases.tsv";
		final String fasta=root+"/members.fa", assignments=root+"/assignments.tsv", references=root+"/models.fa";
		final String[] inputs={"A", "D"};
		final String[] ids={"RF00174_input_A", "RF00174_input_D"};
		final String[] names={"cobalamin_RF00174_bonus_zero", "cobalamin_RF00174_bonus_kept"};
		final byte[][] dna={Prok16sPartialRescueTest.random(160, 101), Prok16sPartialRescueTest.random(160, 202)};
		try(Output m=new Output(manifest); Output p=new Output(panel); Output s=new Output(sources);
			Output a=new Output(aliases); Output f=new Output(fasta); Output assignment=new Output(assignments);
			Output ref=new Output(references)){
			m.writer.println(NcrnaDdlHoldout.MANIFEST_HEADER);
			p.writer.println(NcrnaDdlHoldout.PANEL_HEADER);
			s.writer.println(NcrnaDdlHoldout.SOURCE_HEADER);
			a.writer.println(NcrnaConstructionMembers.GROUP_HEADER);
			assignment.writer.println("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus");
			for(int i=0; i<inputs.length; i++){
				final String id=inputs[i];
				m.writer.println("assembly_"+id+"\t"+id+"\tFIXTURE\t"+root+"/"+id+".fa\t"+root+"/"+id
					+".gff\t.\tACC_"+id+"\t"+root+"/source_"+id+".fa\t.\tFASTA\tNCRNA\tCORE");
				a.writer.println("FIXTURE\t"+id+"\t"+root+"/source_"+id+".fa\t"+root+"/"+id+".fa\tACC_"+id
					+"\t"+(i+1)+"\t0\tfalse\tfalse\tRESOLVED\t"+(i+1)+"\t0\tfalse\tRESOLVED");
				s.writer.println(ids[i]+"\tRF00174\tassembly_"+id+"\t"+id+"\tcontig\t1\t160\t+\t160\t.\t0");
				assignment.writer.println(ids[i]+" details=kept\t0\tRF00174\tRF00174_"+i+"\t"+names[i]+"\tkept");
				final String sequence=new String(dna[i], StandardCharsets.US_ASCII);
				f.writer.println(">"+ids[i]+" details=kept\n"+sequence);
				ref.writer.println(">"+names[i]+"\n"+sequence);
			}
			p.writer.println("0\tassembly_A\tA\t"+root+"/A.fa\t"+root+"/A.gff");
		}
		final boolean amino=Shared.AMINO_IN, iupac=Read.IUPAC_TO_N, uracil=Read.U_TO_T;
		final boolean dots=Read.DOT_DASH_X_TO_N, lower=Read.LOWER_CASE_TO_N;
		try{
			Shared.AMINO_IN=false;
			Read.IUPAC_TO_N=false;
			Read.U_TO_T=false;
			Read.DOT_DASH_X_TO_N=false;
			Read.LOWER_CASE_TO_N=false;
			for(final int bins:new int[]{512, 1024}){
				final int bonus=bins==512 ? 60 : 100;
				final String plain=root+"/plain_"+bins+".tsv", explicit=root+"/explicit_"+bins+".tsv";
				buildBonusFixture(root, plain, bins, new String[0]);
				buildBonusFixture(root, explicit, bins, new String[]{"consensusbonus=0", "consensusrule=none"});
				check(sameBytes(plain, explicit), "Omitted and explicit no-bonus CLI settings must produce identical legacy member-only bytes");
				final NcrnaDdlProfileResource baseline=NcrnaDdlProfileResource.load(plain, names, dna);
				final float[] before=new float[2];
				baseline.newScorer().score(dna[0], before);
				check(baseline.members(0)==0 && baseline.members(1)==1 && before[0]==0,
					"The actual CLI control excludes the held-out voter and leaves its member-only model empty");
				for(final String rule:new String[]{"maxima", "allkeys"}){
					final String path=root+"/"+rule+"_"+bins+".tsv";
					buildBonusFixture(root, path, bins, new String[]{"consensusbonus="+bonus, "consensusrule="+rule});
					final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, dna);
					final float[] after=new float[2];
					resource.newScorer().score(dna[0], after);
					check(resource.members(0)==0 && resource.members(1)==1 && after[0]>0
						&& resource.orderingBonus()==bonus && resource.orderingBonusRule()==NcrnaDdlProfileResource.parseBonusRule(rule),
						"Nondefault CLI bonus reaches the loaded native query index without changing real holdout/member counts");
				}
				reject(()->buildBonusFixture(root, root+"/invalid_"+bins+".tsv", bins,
					new String[]{"consensusbonus="+bonus, "consensusrule=none"}),
					"A positive bonus with no declared rule must fail before profile output");
			}
		}finally{
			Shared.AMINO_IN=amino;
			Read.IUPAC_TO_N=iupac;
			Read.U_TO_T=uracil;
			Read.DOT_DASH_X_TO_N=dots;
			Read.LOWER_CASE_TO_N=lower;
		}
	}

	private static void buildBonusFixture(final String root, final String out, final int bins, final String[] extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(
			"orderonly=t", "in="+root+"/members.fa", "assignments="+root+"/assignments.tsv", "ref="+root+"/models.fa",
			"manifest="+root+"/manifest.tsv", "panel="+root+"/panel.tsv", "sources="+root+"/sources.tsv",
			"aliases="+root+"/aliases.tsv", "expectedinput=2", "expectedmodels=2", "expectedmembers=1",
			"expectedpanel=1", "trainsample=1", "partitionbins="+bins, "paddingblacklist=none", "seed=12345",
			"out="+out, "t=1", "ow=f"));
		args.addAll(Arrays.asList(extra));
		BuildNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));
	}

	private static boolean sameBytes(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] line; (line=a.reader.nextLine())!=null;){
				if(!Arrays.equals(line, b.reader.nextLine())){return false;}
			}
			return b.reader.nextLine()==null;
		}
	}

	/** Three originals produce five vote edges before holdout: two cluster edges
	 * and three singleton edges. The held-out original loses both of its edges. */
	private static void minor(final String root, final String manifest, final String panel, final String aliases){
		final String[] ids={"RF01482_input_P", "RF01482_input_K", "RF01482_input_O"};
		final String[] dna={"ACGTACGT", "GGGGAACC", "TTTTCCAA"};
		final String[] names={"cobalamin_RF01482_c0", "cobalamin_RF01482_extra_"+ids[0],
			"cobalamin_RF01482_extra_"+ids[1], "cobalamin_RF01482_extra_"+ids[2]};
		final byte[][] references={bytes(dna[1]), bytes(dna[0]), bytes(dna[1]), bytes(dna[2])};
		final String sources=root+"/minor.sources.tsv", assignments=root+"/minor.assignments.tsv";
		final String fasta=root+"/minor.fa", additions=root+"/minor.additions.tsv";
		try(Output source=new Output(sources); Output assignment=new Output(assignments);
				Output sequence=new Output(fasta); Output addition=new Output(additions)){
			source.writer.println(NcrnaDdlHoldout.SOURCE_HEADER);
			assignment.writer.println("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus");
			addition.writer.println("locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction");
			for(int i=0; i<ids.length; i++){
				final char genome=i==0 ? 'A' : i==1 ? 'D' : 'C';
				source.writer.println(ids[i]+"\tRF01482\tassembly_"+genome+"\t"+genome+"\tcontig\t1\t8\t+\t8\t.\t0");
				assignment.writer.println(ids[i]+" tid=0\t0\tRF01482\tRF01482_"+i+"\t"
					+(i<2 ? names[0]+"\tkept" : "\torphan"));
				sequence.writer.println(">"+ids[i]+" tid=0\n"+dna[i]);
				addition.writer.println(ids[i]+"\t"+names[i+1]+"\tRF01482\tRF01482\t8\tEXACT_SINGLE_MEMBER_NO_REASSIGNMENT");
			}
			// A foreign member is ignored by this member's source binding.
			source.writer.println("RF00174_input_foreign\tRF00174\tassembly_A\tA\tcontig\t1\t8\t+\t8\t.\t0");
		}
		final NcrnaDdlHoldout holdout=new NcrnaDdlHoldout(manifest, panel, sources, aliases, "RF01482");
		final NcrnaDdlMembers members=NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, additions, names, references, 3, holdout);
		check(members.inputMembers()==3 && members.keptMembers()==1 && members.excludedMembers()==1
			&& members.orphanMembers()==1 && members.addedMembers()==3 && members.excludedAddedMembers()==1
			&& members.votingMembers()==3,
			"Unique original dispositions and model-vote edges must retain their distinct literal conservation totals");
		check(members.members(0).length==1 && members.members(1).length==0
			&& members.members(2).length==1 && members.members(3).length==1
			&& members.members(0)[0]==members.members(2)[0],
			"One source keeps both its cluster and singleton votes using the same read-only bases; excluded singleton stays empty");
		for(final String fault:new String[]{"repeated", "unknown", "length", "sequence"}){
			final String bad=root+"/minor_"+fault+".additions.tsv";
			try(Output out=new Output(bad)){
				out.writer.println("locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction");
				for(int i=0; i<ids.length; i++){
					final String source=fault.equals("unknown") && i==0 ? "RF01482_input_absent" : ids[i];
					final int length=fault.equals("length") && i==0 ? 9 : 8;
					final String row=source+"\tcobalamin_RF01482_extra_"+source+"\tRF01482\tRF01482\t"+length+"\tEXACT_SINGLE_MEMBER_NO_REASSIGNMENT";
					out.writer.println(row);
					if(fault.equals("repeated") && i==0){out.writer.println(row);}
				}
			}
			final byte[][] bound=references.clone();
			if(fault.equals("sequence")){bound[1]=bytes("ACGTACGA");}
			reject(()->NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, bad, names, bound, 3,
				new NcrnaDdlHoldout(manifest, panel, sources, aliases, "RF01482")),
				"Malformed singleton "+fault+" must fail instead of introducing an invented or repeated vote");
		}
		holdout.writeAudit(root+"/minor.audit.tsv");
		try(Input in=new Input(root+"/minor.audit.tsv")){
			in.reader.nextLine();
			final parse.LineParser1 p=new parse.LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.parseInt(2)==3 && p.parseInt(5)==1 && p.parseInt(6)==0,
				"Observe originals once and filter a reused source once; the audit counts excluded unique voters, not duplicate edges");
		}
		reject(()->NcrnaDdlHoldoutTest.loadForeignMinor(fasta, assignments, additions, names, references,
			manifest, panel, sources, aliases), "A mismatched holdout member must fail instead of omitting all minor observations");
	}

	private static void loadForeignMinor(final String fasta, final String assignments, final String additions,
			final String[] names, final byte[][] references, final String manifest, final String panel,
			final String sources, final String aliases){
		NcrnaDdlMembers.loadSubsetRoster(fasta, assignments, additions, names, references, 3,
			new NcrnaDdlHoldout(manifest, panel, sources, aliases));
	}


	/** Historical voters retain real IDs, full lengths, alias and duplicate exclusion. */
	private static void protectedProfiles(final String root, final String manifest, final String panel, final String aliases){
		final String dir=root+"/protected_profiles";
		check(new File(dir).mkdir(), "Protected profile fixture must be fresh");
		final String old="cobalamin_RF01689_old", added="cobalamin_RF01689_added";
		final String[] ids={"RF01689_input_A_locus_1", "RF01689_input_D_locus_2", "RF01689_legacy_alias", "RF01689_legacy_duplicate", "RF01689_legacy_kept"};
		final String[] dna={"ACGTACGTACGTACGTACGTACGT", "CAGTCAGTCAGTCAGTCAGTCAGT", "GGATGGATGGATGGATGGATGGAT", "ACGTACGTACGTACGTACGTACGT", "TGCCTGCCTGCCTGCCTGCCTGCC"};
		final String[] assembly={"A", "D", "B", "C", "D"};
		final String assignmentHeader="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n";
		final String priorHeader="reference_id\tinput_id\tassembly_id\treference_member\tcontig\tstrand\tcore_start1\tcore_stop1\tsource_start1\tsource_stop1\tcontig_length\toriginal_id\toriginal_header\toriginal_taxid\tprior_line\tdecision_line\tsource_kind\tsource_authority\trow_authority\toriginal_genome\tsource_fasta\n";
		final StringBuilder fasta=new StringBuilder(), assignments=new StringBuilder(assignmentHeader), sources=new StringBuilder(NcrnaDdlHoldout.SOURCE_HEADER+'\n'), priors=new StringBuilder(priorHeader);
		for(int i=0; i<ids.length; i++){
			final String header=ids[i]+" details=kept";
			fasta.append('>').append(header).append('\n').append(dna[i]).append('\n');
			assignments.append(header).append("\t-1\tunknown\tRF01689_group\t").append(i==1 ? old : added).append("\tkept\n");
			if(i<2){
				sources.append(ids[i]).append("\tRF01689\tassembly_").append(assembly[i]).append('\t').append(assembly[i])
					.append("\tctg\t1\t24\t+\t24\t.\t0\n");
			}else{
				priors.append("ref_").append(i).append('\t').append(assembly[i]).append("\tassembly_").append(assembly[i])
					.append("\tRF01689\tctg\t+\t2\t23\t1\t24\t100\t").append(ids[i]).append('\t').append(header)
					.append("\t.\t2\t2\tPROTECTED_PRIOR\tfixture\tNA\thistorical_").append(assembly[i]).append('\t').append(root).append('/').append(assembly[i]).append(".fa\n");
			}
		}
		writeText(dir+"/sources.tsv", sources.toString());
		writeText(dir+"/prior.tsv", priors.toString());
		writeText(dir+"/members.fa", fasta.toString());
		writeText(dir+"/assignments.tsv", assignments.toString());
		writeText(dir+"/models.fa", '>'+added+'\n'+dna[4]+'\n');
		final NcrnaDdlHoldout filter=new NcrnaDdlHoldout(manifest, panel, dir+"/sources.tsv", aliases, "RF01689", dir+"/prior.tsv", 3);
		check(filter.boundHistoricalSource(ids[4]+" details=kept") && !filter.boundHistoricalSource(ids[4]+" changed")
			&& !filter.boundHistoricalSource("RF01689_unknown"), "Only exact declared historical headers are admitted");
		final NcrnaDdlMembers members=NcrnaDdlMembers.loadSubsetRoster(dir+"/members.fa", dir+"/assignments.tsv",
			new String[]{added}, new byte[][]{bytes(dna[4])}, 5, filter);
		check(members.inputMembers()==5 && members.assignedMembers(0)==4 && members.keptMembers()==1
			&& members.excludedMembers()==3 && members.removedMembers()==1 && Arrays.equals(members.members(0)[0], bytes(dna[4])),
			"Native panel, historical alias and historical exact duplicate are excluded; one historical voter survives and native removed source is still observed");
		final String[] args={"orderonly=t", "in="+dir+"/members.fa", "assignments="+dir+"/assignments.tsv", "ref="+dir+"/models.fa",
			"manifest="+manifest, "panel="+panel, "sources="+dir+"/sources.tsv", "aliases="+aliases,
			"expectedinput=5", "expectedmodels=1", "expectedmembers=1", "expectedpanel=1", "partitionbins=1024", "paddingblacklist=none",
			"seed=12345", "out="+dir+"/profiles.tsv", "protectedpriors="+dir+"/prior.tsv", "expectedpriors=3", "trd=f", "t=1"};
		BuildNcrnaDdlProfiles.main(args);
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(dir+"/profiles.tsv", new String[]{added}, new byte[][]{bytes(dna[4])});
		check(resource.members(0)==1, "Actual ordering-only CLI retains the historical voter in its native profile format");
		writeText(dir+"/native.fa", '>'+ids[0]+" details=kept\n"+dna[0]+"\n>"+ids[1]+" details=kept\n"+dna[1]+'\n');
		writeText(dir+"/native.assignments.tsv", assignmentHeader+ids[0]+" details=kept\t-1\tunknown\tRF01689_old\t"+old+"\tkept\n"
			+ids[1]+" details=kept\t-1\tunknown\tRF01689_old\t"+old+"\tkept\n");
		writeText(dir+"/old.models.fa", '>'+old+'\n'+dna[1]+'\n');
		final String[] oldArgs=Arrays.copyOf(args, 16);
		oldArgs[1]="in="+dir+"/native.fa";
		oldArgs[2]="assignments="+dir+"/native.assignments.tsv";
		oldArgs[3]="ref="+dir+"/old.models.fa";
		oldArgs[8]="expectedinput=2";
		oldArgs[15]="out="+dir+"/old.profiles.tsv";
		BuildNcrnaDdlProfiles.main(oldArgs);
		BuildNcrnaDdlProfiles.main(new String[]{"appendprofiles=t", "in="+dir+"/old.profiles.tsv", "in2="+dir+"/profiles.tsv",
			"ref="+dir+"/old.models.fa", "ref2="+dir+"/models.fa", "expectedmodels=1,1", "out="+dir+"/appended.tsv"});
		final NcrnaDdlProfileResource oldResource=NcrnaDdlProfileResource.load(dir+"/old.profiles.tsv", new String[]{old}, new byte[][]{bytes(dna[1])});
		final NcrnaDdlProfileResource appended=NcrnaDdlProfileResource.load(dir+"/appended.tsv", new String[]{old, added}, new byte[][]{bytes(dna[1]), bytes(dna[4])});
		check(appended.members(0)==1 && appended.members(1)==1, "Protected construction rows append to ordinary profiles without changing voter counts");
		for(final String query:dna){
			final float[] oldScore=new float[1], newScore=new float[1], joined=new float[2];
			oldResource.newScorer().score(bytes(query), oldScore);
			resource.newScorer().score(bytes(query), newScore);
			appended.newScorer().score(bytes(query), joined);
			check(Float.floatToIntBits(joined[0])==Float.floatToIntBits(oldScore[0]) && Float.floatToIntBits(joined[1])==Float.floatToIntBits(newScore[0]),
				"Append preserves both old and protected-prior profile query scores exactly");
		}

		for(int fault=0; fault<6; fault++){
			final String prefix=dir+"/bad"+fault;
			final String[] bad=args.clone();
			bad[15]="out="+prefix+".tsv";
			if(fault==0){
				writeText(prefix+".prior.tsv", "malformed_companion_header\n");
				bad[16]="protectedpriors="+prefix+".prior.tsv";
			}else if(fault==1){bad[17]="expectedpriors=2";}
			else if(fault==2){
				writeText(prefix+".prior.tsv", priors.toString().replace("\tassembly_B\t", "\tassembly_A\t"));
				bad[16]="protectedpriors="+prefix+".prior.tsv";
			}else if(fault==3){
				writeText(prefix+".prior.tsv", priors.toString().replace(root+"/D.fa", root+"/wrong.fa"));
				bad[16]="protectedpriors="+prefix+".prior.tsv";
			}else{
				writeText(prefix+".fa", fasta.toString().replace(fault==4 ? ids[4]+" details=kept" : dna[4], fault==4 ? ids[4]+" changed" : dna[4]+"A"));
				bad[1]="in="+prefix+".fa";
			}
			reject(()->BuildNcrnaDdlProfiles.main(bad), "Malformed historical source binding must fail before writing a profile");
			check(!new File(prefix+".tsv").exists(), "Invalid historical input cannot publish profile rows");
		}
		reject(()->NcrnaDdlMembers.loadSubsetRoster(dir+"/members.fa", dir+"/assignments.tsv", new String[]{added}, new byte[][]{bytes(dna[4])}, 5,
			new NcrnaDdlHoldout(manifest, panel, dir+"/sources.tsv", aliases, "RF01689")),
			"The default native-only route still rejects unbound historical IDs");
	}
	private static void writeText(final String path, final String contents){try(Output out=new Output(path)){out.writer.print(contents);}}

	private static void writeFasta(final String path, final String[] ids, final String[] dna){
		try(Output out=new Output(path)){
			// The cross-assembly duplicate precedes the panel gene deliberately.
			for(final int i:new int[]{2, 3, 4, 0, 1}){
				out.writer.println(">"+ids[i]+" details=kept");
				out.writer.println(dna[i]);
			}
		}
	}

	private static void reject(final Runnable action, final String why){
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, why);
	}

	private static void literal(){
		final NcrnaExactSequenceSet exact=new NcrnaExactSequenceSet(true);
		exact.add(bytes("ACGRNT"));
		check(exact.contains(bytes("acgrnu")) && !exact.contains(bytes("ACGANT")),
			"Literal IUPAC equality folds case/U but never treats ambiguous symbols as overlapping alternatives");
		boolean rejected=false;
		try{new NcrnaExactSequenceSet().add(bytes("ACGRNT"));}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "The existing defined-base-only exact-set contract must remain unchanged by the new explicit mode");
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
}
