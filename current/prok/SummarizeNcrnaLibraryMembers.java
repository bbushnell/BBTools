package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.IntObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.DoubleList;
import structures.IntHashMap;
import structures.IntList;

/** Joins saved construction members, assignments and consensus lengths.
 * Uses A_SampleByteFile's parser/ByteFile/output lifecycle and the existing
 * native FASTA stream owner. No source genome, alignment or caller is read.
 * Construction votes and original-coordinate gene multiplicity stay separate.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaLibraryMembers{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaLibraryMembers.class, false);
		try{
			final SummarizeNcrnaLibraryMembers tool=new SummarizeNcrnaLibraryMembers(pp.args);
			tool.process();
			pp.outstream.println("NCRNA_LIBRARY_MEMBER_SHAPES_COMPLETE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaLibraryMembers(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		String names=null, models=null, assignments=null;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("originals")){originals=b;}
			else if(a.equals("labels")){names=b;}
			else if(a.equals("consensus")){models=b;}
			else if(a.equals("assignments")){assignments=b;}
			else if(a.equals("expectedmembers")){expectedMembers=Parse.parseIntKMG(b);}
			else if(a.equals("expectedoriginals")){expectedOriginals=Parse.parseIntKMG(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown library member parameter: "+arg);}
		}
		membersPath=parser.in1;
		prefix=parser.out1;
		require(membersPath!=null && originals!=null && prefix!=null && names!=null && models!=null && assignments!=null
			&& expectedMembers>0 && expectedOriginals>0 && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Use saved member FASTA, originals, matching labels/consensus/assignments lists and exact positive populations; no sampling");
		labels=names.split(",", -1);
		consensusPaths=models.split(",", -1);
		assignmentPaths=assignments.split(",", -1);
		require(labels.length>0 && labels.length==consensusPaths.length && labels.length==assignmentPaths.length,
			"Each saved library needs one label, consensus and assignment file");
		for(int i=0; i<labels.length; i++){
			require(labels[i].matches("[A-Za-z0-9_]+"), "Labels must be distinct plain file-name components");
			for(int j=0; j<i; j++){require(!labels[i].equals(labels[j]), "Library labels must not overwrite each other");}
			RrnaResourceIO.fresh(prefix+"."+labels[i], ".summary.tsv", ".models.tsv");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Metadata         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		Shared.TRIM_READ_DESCRIPTION=false;
		NcrnaDdlSavedData.sequences(membersPath, read->{
			final int id=memberId(read.id);
			final int from=attribute(read.id, " cmfrom=");
			final int to=attribute(read.id, " cmto="), length=attribute(read.id, " cmlen=");
			require(!members.contains(id) && from>=1 && to>=from && length>=to && read.length()>0,
				"Construction representatives need unique IDs, real model coordinates and actual saved sequence lengths");
			final Member member=new Member(read.length(), from, to, length,
				booleanAttribute(read.id, " incomplete_model_coverage="));
			members.put(id, member);
			memberOrder.add(member);
		});
		require(members.size()==expectedMembers, "The member-length denominator must equal the complete construction population");
		loadOriginals();
		for(int i=0; i<labels.length; i++){library(i);}
	}

	private void loadOriginals(){
		final LineParser1 p=new LineParser1('\t');
		final IntHashMap seen=new IntHashMap();
		try(NcrnaTextIO.Input in=new NcrnaTextIO.Input(originals)){
			RrnaResourceIO.header(in.reader, PrepareNcrnaRosterBaseline.ORIGINAL_HEADER, originals);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==29 && p.parseInt(0)>=0 && !seen.contains(p.parseInt(0)) && p.termEquals("true", 7)
					&& p.parseInt(18)>0 && p.parseInt(19)>=p.parseInt(18) && p.parseInt(20)>=p.parseInt(19),
					"Original gene rows must be unique and retain the full saved 29-column schema");
				seen.put(p.parseInt(0), 1);
				final int representative=p.parseInt(14);
				require(representative<0 || members.contains(representative), "Every eligible original gene must bind a saved construction representative");
				require(p.termEquals("true", 21) || p.termEquals("false", 21), "Original coverage flags must retain their saved boolean meaning");
				final Member own=members.get(p.parseInt(13));
				if(own!=null && p.parseInt(13)==representative && own.from==p.parseInt(18)
						&& own.to==p.parseInt(19) && own.modelLength==p.parseInt(20)
						&& own.length==p.parseInt(5)-p.parseInt(4)+1 && own.incomplete==p.termEquals("true", 21)
						&& p.parseInt(17)==0 && p.termEquals("ELIGIBLE_NATIVE_HIT", 22)){own.originalBound=true;}
				if(p.parseInt(18)-1!=44 || p.parseInt(20)-p.parseInt(19)!=20){continue;}
				originalShort++;
				if(representative<0){originalShortExcluded++;}
				else{originalShortRepresentatives.add(representative);}
			}
		}
		require(seen.size()==expectedOriginals && originalShort==originalShortExcluded+originalShortRepresentatives.size,
			"All original genes and every short-class construction disposition must be conserved");
		for(final Member member:memberOrder){
			require(member.originalBound, "Each construction representative needs its own exact original component, native model coordinates and sequence span");
		}
	}

	private void library(final int index){
		final ObjectMap<String, Model> names=new ObjectMap<String, Model>(String.class, Model.class);
		final ArrayList<Model> models=new ArrayList<Model>();
		for(final Member member:memberOrder){member.model=null;}
		NcrnaDdlSavedData.sequences(consensusPaths[index], read->{
			final String name=RrnaResourceIO.first(read.id);
			require(read.length()>0 && !names.contains(name), "Consensus models need unique names and real emitted lengths");
			final Model model=new Model(name, read.length());
			names.put(name, model);
			models.add(model);
		});
		require(!models.isEmpty(), "An accepted saved library must contain consensus models");
		final LineParser1 p=new LineParser1('\t');
		int assigned=0;
		try(NcrnaTextIO.Input in=new NcrnaTextIO.Input(assignmentPaths[index])){
			RrnaResourceIO.header(in.reader, "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus", assignmentPaths[index]);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==6 && p.termStartsWith("RF00174_X_", 0) && p.termEquals("kept", 5),
					"These complete libraries must retain one kept-model assignment for every construction member");
				p.setBounds(0);
				int end=p.a()+10;
				while(end<p.b() && p.line()[end]!=' '){end++;}
				final Member member=members.get(Parse.parseInt(p.line(), p.a()+10, end));
				final Model model=names.get(RrnaResourceIO.first(p.parseString(4)));
				require(member!=null && member.model==null && model!=null, "Each saved member must have exactly one known model assignment");
				member.model=model;
				model.ratios.add((double)member.length/model.length);
				if(member.shortBoth){model.shortMembers++;}
				if(member.full){model.fullMembers++;}
				if(member.incomplete){model.incompleteMembers++;}
				assigned++;
			}
		}
		require(assigned==expectedMembers, "Every library must conserve all construction assignments before its mixture or length spread is reported");
		for(int i=0; i<originalShortRepresentatives.size; i++){
			final Model model=members.get(originalShortRepresentatives.get(i)).model;
			require(model!=null, "An eligible original short gene cannot lose its assigned construction model");
			model.originalShort++;
		}
		write(index, models);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Reports and Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final int index, final ArrayList<Model> models){
		final double[][] quantiles=new double[3][models.size()];
		int shortMembers=0, fullMembers=0, incomplete=0, shortModels=0, pure=0, over90=0;
		int mixedFullModels=0, shortWithFull=0, shortWithOther=0, originalModels=0, originalWithFull=0;
		final ByteBuilder row=new ByteBuilder();
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(prefix+"."+labels[index]+".models.tsv")){
			out.writer.println("model\tconsensus_bases\tmembers\tboth_ends_short_members\texact_full_span_members\tother_members\tincomplete_coverage_members_overlaps_short\toriginal_short_genes\tmember_length_ratio_5th\tmember_length_ratio_median\tmember_length_ratio_95th");
			for(int i=0; i<models.size(); i++){
				final Model model=models.get(i);
				final int n=model.ratios.size;
				require(n>0 && model.shortMembers+model.fullMembers<=n, "Every emitted consensus must retain members and disjoint short/full classes");
				Arrays.sort(model.ratios.array, 0, n);
				for(int q=0; q<3; q++){quantiles[q][i]=quantile(model.ratios.array, n, QUANTILES[q]);}
				shortMembers+=model.shortMembers;
				fullMembers+=model.fullMembers;
				incomplete+=model.incompleteMembers;
				if(model.shortMembers>0){
					shortModels++;
					if(model.shortMembers==n){pure++;}
					if((long)model.shortMembers*10>(long)n*9){over90++;}
					if(model.fullMembers>0){
						mixedFullModels++;
						shortWithFull+=model.shortMembers;
					}
					if(model.shortMembers<n){shortWithOther+=model.shortMembers;}
				}
				if(model.originalShort>0){originalModels++;}
				if(model.fullMembers>0){originalWithFull+=model.originalShort;}
				row.clear().append(model.name).tab().append(model.length).tab().append(n).tab()
					.append(model.shortMembers).tab().append(model.fullMembers).tab().append(n-model.shortMembers-model.fullMembers).tab()
					.append(model.incompleteMembers).tab()
					.append(model.originalShort);
				for(int q=0; q<3; q++){row.tab().append(quantiles[q][i], 8);}
				out.writer.print(row.nl());
			}
		}
		for(final double[] values:quantiles){Arrays.sort(values);}
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(prefix+"."+labels[index]+".summary.tsv")){
			out.writer.println("metric\tvalue");
			metric(out, "models", models.size());
			metric(out, "construction_members", expectedMembers);
			metric(out, "incomplete_coverage_members", incomplete);
			metric(out, "exact_full_span_members", fullMembers);
			metric(out, "both_ends_short_members", shortMembers);
			metric(out, "other_members", expectedMembers-shortMembers-fullMembers);
			metric(out, "models_with_short_members", shortModels);
			metric(out, "short_only_models", pure);
			metric(out, "models_over_90_percent_short", over90);
			metric(out, "short_models_with_full_span_members", mixedFullModels);
			metric(out, "short_members_mixed_with_full_span", shortWithFull);
			metric(out, "short_members_mixed_with_any_other_class", shortWithOther);
			metric(out, "original_genes", expectedOriginals);
			metric(out, "original_both_ends_short_genes", originalShort);
			metric(out, "original_short_genes_outside_construction", originalShortExcluded);
			metric(out, "models_with_original_short_genes", originalModels);
			metric(out, "original_short_genes_in_models_with_full_span_members", originalWithFull);
			for(int q=0; q<3; q++){
				metric(out, "median_across_models_of_member_length_ratio_"+QUANTILES[q]+"th", quantile(quantiles[q], models.size(), 50));
			}
			metric(out, "5th_percentile_across_model_median_length_ratios", quantile(quantiles[1], models.size(), 5));
			metric(out, "95th_percentile_across_model_median_length_ratios", quantile(quantiles[1], models.size(), 95));
		}
	}

	private static double quantile(final double[] sorted, final int size, final int percentile){
		assert(size>0 && size<=sorted.length && percentile>0 && percentile<=100)
			: "Nearest-rank member/model quantiles require a complete nonempty sorted population";
		return sorted[(int)(((long)size*percentile+99)/100)-1];
	}

	private static int memberId(final String id){
		require(id.startsWith("RF00174_X_"), "Use saved RF00174 construction member headers");
		final int end=id.indexOf(' ');
		return Parse.parseInt(id, 10, end<0 ? id.length() : end);
	}

	private static int attribute(final String header, final String name){
		final int from=attributeStart(header, name);
		final int end=header.indexOf(' ', from);
		return Parse.parseInt(header, from, end<0 ? header.length() : end);
	}

	private static boolean booleanAttribute(final String header, final String name){
		final int from=attributeStart(header, name);
		final int space=header.indexOf(' ', from), end=space<0 ? header.length() : space;
		if(end-from==4 && header.startsWith("true", from)){return true;}
		if(end-from==5 && header.startsWith("false", from)){return false;}
		throw new IllegalArgumentException("Saved member coverage must be an explicit true/false value");
	}

	private static int attributeStart(final String header, final String name){
		final int start=header.indexOf(name);
		if(start<0){throw new IllegalArgumentException("Missing saved member metadata: "+name);}
		return start+name.length();
	}

	private void metric(final NcrnaTextIO.Output out, final String key, final double value){
		out.writer.print(summaryRow.clear().append(key).tab().append(value, 8).nl());
	}

	private void metric(final NcrnaTextIO.Output out, final String key, final long value){
		out.writer.print(summaryRow.clear().append(key).tab().append(value).nl());
	}

	private static void require(final boolean condition, final String reason){
		if(!condition){throw new IllegalArgumentException(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Metadata      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Member{
		Member(final int length_, final int from_, final int to_, final int modelLength_, final boolean incomplete_){
			assert(length_>0 && from_>0 && to_>=from_ && modelLength_>=to_)
				: "Member length/cohort metadata must retain valid saved sequence and native model bounds";
			length=length_;
			from=from_;
			to=to_;
			modelLength=modelLength_;
			shortBoth=from-1==44 && modelLength-to==20;
			full=from==1 && to==modelLength;
			incomplete=incomplete_;
		}
		final int length, from, to, modelLength;
		final boolean shortBoth, full, incomplete;
		boolean originalBound;
		Model model;
	}

	private static final class Model{
		Model(final String name_, final int length_){
			name=name_;
			length=length_;
		}
		final String name;
		final int length;
		final DoubleList ratios=new DoubleList(8);
		int shortMembers, fullMembers, incompleteMembers, originalShort;
	}

	private String membersPath, originals, prefix;
	private String[] labels, consensusPaths, assignmentPaths;
	private int expectedMembers, expectedOriginals, originalShort, originalShortExcluded;
	private final IntObjectMap<Member> members=new IntObjectMap<Member>(Member.class);
	private final ArrayList<Member> memberOrder=new ArrayList<Member>();
	private final IntList originalShortRepresentatives=new IntList();
	private final ByteBuilder summaryRow=new ByteBuilder();
	private static final int[] QUANTILES={5,50,95};
}
