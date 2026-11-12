package prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectIntMap;
import map.ObjectMap;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaSubtypeEvidence.Support;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Inventory of existing reference parents, never a caller or an adjudicator.
 * One thread reads canonical GFFs once; original FASTAs are never opened. Family
 * unions count each parent once. Member support and the two cobalamin reference
 * views are deliberately nonadditive. Catalog source memberships remain strata,
 * not biological negative labels or proof of independent organisms.
 * @author Ganyu
 */
final class NcrnaReferenceInventory{

	static void process(final String census, final String memberships, final String reference,
		final String members, final String witnesses, final String events, final String prefix){
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		assert(catalog!=null && catalog.size()>0) : "Inventory requires the complete original-input catalog";
		new NcrnaReferenceInventory(catalog, census, memberships, reference, members, witnesses, events, prefix).run();
	}

	private NcrnaReferenceInventory(final NcrnaEvidenceCatalog catalog_, final String census_, final String memberships_, final String reference_,
		final String members_, final String witnesses_, final String events_, final String prefix_){
		catalog=catalog_;
		census=census_;
		memberships=memberships_;
		reference=reference_;
		members=members_;
		witnesses=witnesses_;
		events=events_;
		prefix=prefix_;
		final ArrayList<String> inputs=new ArrayList<String>(catalog.size());
		for(int i=0; i<catalog.size(); i++){inputs.add(catalog.get(i).id);}
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){
			final Genome genome=catalog.get(binding.id);
			require(genome!=null && genome.fasta.equals(binding.fasta), "Reference must bind its exact original catalog FASTA: "+binding.id);
			bindings.put(binding.id, binding);
		}
		require(bindings.size()==catalog.size(), "Reference population must cover every catalog input, including zero-target inputs");
		evidence=new NcrnaSubtypeEvidence(members, witnesses, events, inputs);
		memberFamily=new int[evidence.size()];
		for(int m=0; m<evidence.size(); m++){
			final String family=evidence.family(m);
			if(!familyIndex.contains(family)){
				familyIndex.put(family, families.size());
				families.add(family);
			}
			memberFamily[m]=familyIndex.get(family);
		}
		familyTotals=stats(families.size());
		memberTotals=stats(evidence.size());
		recoveryTotals=stats(evidence.size());
		exclusiveTotals=stats(evidence.size());
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			for(int j=0; j<genome.memberships(); j++){
				final String source=genome.membership(j).source;
				if(!sourceIndex.contains(source)){
					sourceIndex.put(source, sources.size());
					sources.add(source);
					sourceTotals.add(stats(families.size()));
				}
			}
		}
		lastSourceInput=new int[sources.size()];
		Arrays.fill(lastSourceInput, -1);
	}

	private void run(){
		final long[][] union=new long[families.size()][4];
		final long[][] member=new long[evidence.size()][6];
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		try(Output inputs=new Output(prefix+SUFFIXES[0]); Output membersOut=new Output(prefix+SUFFIXES[1]);
			Output loci=new Output(prefix+SUFFIXES[2])){
			inputs.writer.println(INPUT_HEADER);
			membersOut.writer.println(MEMBER_HEADER);
			loci.writer.println(LOCUS_HEADER);
			for(int i=0; i<catalog.size(); i++){
				for(final long[] counts:union){Arrays.fill(counts, 0);}
				for(final long[] counts:member){Arrays.fill(counts, 0);}
				final Genome genome=catalog.get(i);
				final Binding binding=bindings.get(genome.id);
				require(binding!=null, "Each original input needs its canonical reference pointer");
				final ArrayList<Row> rows=codec.read(binding.gff, genome.id);
				for(final Row row:rows){
					if(!familyIndex.contains(row.locus.family)){continue;}
					final int f=familyIndex.get(row.locus.family);
					final Support support=evidence.support(row);
					union[f][row.mask ? 1 : 0]++;
					if(support.priorType>=0){union[f][2]++;}
					if(!row.mask && support.multi()){union[f][3]++;}
					countMembers(row, support, f, member);
					writeLocus(row, support, loci);
					locusRows++;
				}
				for(int f=0; f<families.size(); f++){
					final long[] count=union[f];
					familyTotals[f].add(count[0], count[1], count[2], count[3]);
					field(buffer.clear(), genome.id).tab();
					field(buffer, families.get(f)).tab();
					appendInputCounts(buffer, count[0], count[1]);
					buffer.tab().append(count[2]).tab().append(count[3]).nl();
					inputs.writer.print(buffer);
				}
				for(int m=0; m<evidence.size(); m++){
					final long[] count=member[m];
					memberTotals[m].add(count[0], count[1], 0, 0);
					final boolean cobalamin=NcrnaSubtypeEvidence.FAMILY.equals(evidence.family(m));
					if(cobalamin){
						recoveryTotals[m].add(count[2], count[3], 0, 0);
						exclusiveTotals[m].add(count[4], count[5], 0, 0);
					}
					field(buffer.clear(), genome.id).tab();
					memberKey(buffer, m);
					buffer.tab().append(count[0]).tab().append(count[1]).tab().append(count[0]>0);
					if(cobalamin){
						for(int c=2; c<6; c++){buffer.tab().append(count[c]);}
					}else{
						buffer.append("\tNA\tNA\tNA\tNA");
					}
					membersOut.writer.print(buffer.nl());
				}
				for(int j=0; j<genome.memberships(); j++){
					final int s=sourceIndex.get(genome.membership(j).source);
					if(lastSourceInput[s]==i){continue;}
					lastSourceInput[s]=i;
					for(int f=0; f<families.size(); f++){
						final long[] count=union[f];
						sourceTotals.get(s)[f].add(count[0], count[1], count[2], count[3]);
					}
				}
			}
		}
		evidence.finish();
		long expected=0;
		for(final Counts counts:familyTotals){
			counts.validate();
			require(counts.inputs==catalog.size(), "Every family includes all original input IDs, including zeros");
			expected=Math.addExact(expected, Math.addExact(counts.scored, counts.masked));
		}
		require(expected==locusRows, "Union totals must conserve distinct retained parent rows, never summed member support");
		writeSummaries();
		writeProvenance();
	}

	private void countMembers(final Row row, final Support item, final int family, final long[][] counts){
		assert(memberFamily.length==counts.length) : "Member counters follow the explicit support map";
		for(int m=0; m<counts.length; m++){
			if(memberFamily[m]!=family){continue;}
			final int bit=1<<m;
			if((item.scored&bit)!=0){counts[m][0]++;}
			if((item.marginal&bit)!=0){counts[m][1]++;}
			if(!NcrnaSubtypeEvidence.FAMILY.equals(row.locus.family)
				|| !NcrnaSubtypeEvidence.includedInView(row, item, m)){continue;}
			counts[m][NcrnaSubtypeEvidence.maskedInView(row, item, m, false) ? 3 : 2]++;
			counts[m][NcrnaSubtypeEvidence.maskedInView(row, item, m, true) ? 5 : 4]++;
		}
	}

	private void writeLocus(final Row row, final Support item, final Output output){
		assert(row.input.equals(item.input)) : "Support must retain the original parent binding";
		field(buffer.clear(), row.input).tab();
		field(buffer, row.id).tab();
		field(buffer, row.locus.family).tab();
		field(buffer, row.locus.seqid).tab().append(row.locus.start).tab().append(row.locus.stop).tab().append(row.locus.strand).tab();
		field(buffer, row.mask ? "dont_care" : "scored").tab();
		field(buffer, row.reason==null ? "." : row.reason).tab();
		field(buffer, row.model==null ? "." : row.model).tab();
		appendMembers(buffer, item.scored);
		buffer.tab();
		appendMembers(buffer, item.marginal);
		buffer.tab().append(item.priorType>=0).tab().append(item.multi()).tab().append(item.bridgeWitnesses()).tab();
		field(buffer, row.evidence).nl();
		output.writer.print(buffer);
	}

	private void appendMembers(final ByteBuilder out, final int bits){
		assert(bits>=0) : "Member support uses positive bits only";
		if(bits==0){
			out.append('.');
			return;
		}
		boolean first=true;
		for(int m=0; m<evidence.size(); m++){
			if((bits&(1<<m))==0){continue;}
			if(!first){out.append(',');}
			field(out, evidence.name(m));
			first=false;
		}
	}

	private void writeSummaries(){
		try(Output out=new Output(prefix+SUFFIXES[3]); Output strata=new Output(prefix+SUFFIXES[4])){
			out.writer.println(SUMMARY_HEADER);
			strata.writer.println("source\t"+SUMMARY_HEADER);
			for(int f=0; f<families.size(); f++){
				summaryRow(buffer.clear(), "UNION", families.get(f), ".", familyTotals[f]);
				out.writer.print(buffer);
			}
			for(int m=0; m<evidence.size(); m++){
				summaryRow(buffer.clear(), "QUALIFIED_MEMBER_SUPPORT", evidence.family(m), evidence.name(m), memberTotals[m]);
				out.writer.print(buffer);
				if(NcrnaSubtypeEvidence.FAMILY.equals(evidence.family(m))){
					summaryRow(buffer.clear(), "RECOVERY", evidence.family(m), evidence.name(m), recoveryTotals[m]);
					out.writer.print(buffer);
					summaryRow(buffer.clear(), "EXCLUSIVE", evidence.family(m), evidence.name(m), exclusiveTotals[m]);
					out.writer.print(buffer);
				}
			}
			for(int s=0; s<sources.size(); s++){
				for(int f=0; f<families.size(); f++){
					field(buffer.clear(), sources.get(s)).tab();
					summaryRow(buffer, "UNION", families.get(f), ".", sourceTotals.get(s)[f]);
					strata.writer.print(buffer);
				}
			}
		}
	}

	private void summaryRow(final ByteBuilder out, final String view, final String family, final String member, final Counts c){
		c.validate();
		field(out, view).tab();
		field(out, family).tab();
		field(out, member).tab();
		out.append(c.inputs).tab().append(c.positive).tab().append(c.maskOnly).tab().append(c.zero)
			.tab().append(c.scored).tab().append(c.masked).tab().append(c.protectedParents).tab().append(c.multi)
			.append("\tORIGINAL_INPUT_IDS_NOT_PROVEN_INDEPENDENT_ORGANISMS;MEMBER_VIEWS_NONADDITIVE\n");
	}

	/** Written last, only after all payload writers and conservation checks succeed. */
	private void writeProvenance(){
		try(Output out=new Output(prefix+SUFFIXES[5])){
			out.writer.println("role\tpath_or_value");
			for(final String[] pair:new String[][]{{"census", census}, {"memberships", memberships}, {"reference", reference}, {"members", members}, {"witnesses", witnesses},
				{"events", events}, {"status", "REFERENCE_INVENTORY_COMPLETE_NO_SELECTION_OR_ADJUDICATION"}}){
				field(buffer.clear(), pair[0]).tab();
				field(buffer, pair[1]).nl();
				out.writer.print(buffer);
			}
		}
	}

	private void memberKey(final ByteBuilder out, final int m){
		assert(m>=0 && m<evidence.size()) : "Member metadata indexes the same validated support map";
		field(out, evidence.family(m)).tab();
		field(out, evidence.name(m)).tab();
		field(out, evidence.referenceModel(m));
	}
	private static void appendInputCounts(final ByteBuilder out, final long scored, final long masked){
		assert(scored>=0 && masked>=0) : "Reference counts cannot be negative";
		out.append(scored).tab().append(masked).tab().append(scored>0).tab().append(scored==0 && masked>0).tab().append(scored==0 && masked==0);
	}
	private static Counts[] stats(final int length){
		assert(length>0) : "A declared family/member population must not be empty";
		final Counts[] array=new Counts[length];
		for(int i=0; i<length; i++){array[i]=new Counts();}
		return array;
	}
	private static ByteBuilder field(final ByteBuilder out, final String value){return BuildNcrnaReference.field(out, value);}
	private static void require(final boolean ok, final String why){
		if(!ok){throw new IllegalArgumentException(why);}
	}
	private static final class Counts{
		void add(final long s, final long m, final long p, final long multi_){
			assert(s>=0 && m>=0 && p>=0 && multi_>=0 && p<=s && multi_<=s) : "Distinct scored parents bound protected and multi-supported counts";
			inputs++;
			if(s>0){
				positive++;
			}else if(m>0){
				maskOnly++;
			}else{
				zero++;
			}
			scored=Math.addExact(scored, s);
			masked=Math.addExact(masked, m);
			protectedParents=Math.addExact(protectedParents, p);
			multi=Math.addExact(multi, multi_);
		}
		void validate(){require(inputs==positive+maskOnly+zero, "Positive, masked-only and zero-target inputs must partition the input population");}
		long inputs, positive, maskOnly, zero, scored, masked, protectedParents, multi;
	}

	private final NcrnaEvidenceCatalog catalog;
	private final NcrnaSubtypeEvidence evidence;
	private final String census, memberships, reference, members, witnesses, events, prefix;
	private final ObjectMap<String, Binding> bindings=new ObjectMap<String, Binding>(String.class, Binding.class);
	private final ObjectIntMap<String> familyIndex=new ObjectIntMap<String>(String.class), sourceIndex=new ObjectIntMap<String>(String.class);
	private final ArrayList<String> families=new ArrayList<String>(), sources=new ArrayList<String>();
	private final ArrayList<Counts[]> sourceTotals=new ArrayList<Counts[]>();
	private final Counts[] familyTotals, memberTotals, recoveryTotals, exclusiveTotals;
	private final int[] memberFamily, lastSourceInput;
	private final ByteBuilder buffer=new ByteBuilder();
	private long locusRows;
	static final String[] SUFFIXES={".inputs.tsv", ".members.tsv", ".loci.tsv.gz", ".summary.tsv", ".strata.tsv", ".provenance.tsv"};
	static final String INPUT_HEADER="input_id\tfamily\tscored_parents\tmasked_parents\tpositive_input\tmasked_only_input\tzero_target_input\tprotected_parents\tmulti_supported_scored_parents";
	static final String MEMBER_HEADER="input_id\tfamily\tmember\treference_model\tqualified_scored_parents\tqualified_marginal_parents\tpositive_input\trecovery_scored_parents\trecovery_masked_parents\texclusive_scored_parents\texclusive_masked_parents";
	static final String LOCUS_HEADER="input_id\tparent_id\tfamily\tseqid\tstart\tstop\tstrand\tstatus\tmask_reason\trepresentative_model\tscored_members\tmarginal_members\tprotected_prior\tmulti_supported\tbridge_witnesses\tevidence";
	static final String SUMMARY_HEADER="view\tfamily\tmember\tinputs\tpositive_inputs\tmasked_only_inputs\tzero_target_inputs\tscored_parents\tmasked_parents\tprotected_parents\tmulti_supported_scored_parents\tqualification";
}
