package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Summarizes original winners of already classified absent-call losses.
 * Adapted from template.A_SampleByteFile, with native byte parsing and output.
 * Each input is a completed CompareCobalaminBenchmarkGrades locus table; panel
 * groups remain separate and the UNION model counts deduplicate across panels.
 * Locus counts in UNION are evidence rows, not deduplicated biological features.
 * Does not repeat matching, classify overlaps or export/restore any library.
 * Full model text survives beside native first-token IDs. Already-retained
 * winners are listed separately because adding them again cannot fix absence.
 * Instances are single-use and single-threaded; retained objects are per model.
 * @author Raiden
 */
public final class SummarizeNcrnaLostModels{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaLostModels.class, false);
		try{new SummarizeNcrnaLostModels(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaLostModels(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		String panels=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("panels")){panels=b;}
			else{require(parser.parse(arg, a, b), "Unknown lost-model summary parameter: "+arg);}
		}
		require(parser.in1!=null && panels!=null && parser.out1!=null && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=completed_locus_tables panels=matching_labels out=fresh_prefix without sampling or overwrite");
		inputs=parser.in1.split(",");
		labels=panels.split(",");
		out=parser.out1;
		require(inputs.length==labels.length && inputs.length>0, "Each completed input table needs one explicit panel label");
		for(int i=0; i<labels.length; i++){
			require(labels[i].matches("[A-Za-z0-9_-]+") && !labels[i].equals("UNION"), "Use a nonempty plain panel label distinct from UNION");
			for(int j=0; j<i; j++){
				require(!labels[i].equals(labels[j]) && !inputs[i].equals(inputs[j]), "Duplicate panel label or input would double-count evidence: "+labels[i]);
			}
		}
		counts=new long[inputs.length+1][4][6];
		require(Tools.testInputFiles(false, true, inputs) && Tools.testOutputFiles(false, false, false,
			out+".models.tsv", out+".removed.ids.txt", out+".already_present.loci.tsv", out+".summary.tsv",
			out+".RF00174.ids.txt", out+".RF01482.ids.txt", out+".RF01689.ids.txt"),
			"Input tables must be readable and all summaries fresh");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads each completed classification once; does not inspect DNA or call GFFs. */
	private void process(){
		try(Output present=new Output(out+".already_present.loci.tsv")){
			present.writer.println("panel\t"+HEADER);
			for(int panel=0; panel<inputs.length; panel++){read(panel, present);}
		}
		write();
	}

	/** Retains model identities only for absence-class losses. */
	private void read(final int panel, final Output presentOut){
		assert(panel>=0 && panel<inputs.length) : "Panel index selects the matching input, label and evidence counter axis";
		try(Input in=new Input(inputs[panel])){
			RrnaResourceIO.header(in.reader, HEADER, inputs[panel]);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==13, "Completed locus classification requires thirteen fields");
				if(p.termEquals("GAINED_GRADED_REFERENCE", 6)){continue;}
				require(p.termEquals("LOST_GRADED_REFERENCE", 6), "Expected a recorded graded loss or gain");
				final int member=member();
				increment(panel, member, 0);
				if(p.termEquals("OVERLAPPING_CLAN_CALL_PRESENT", 7)){continue;}
				require(p.termEquals("NO_OVERLAPPING_CLAN_CALL", 7), "Unknown loss-presence class");
				increment(panel, member, 1);
				final boolean present=p.termEquals("true", 10);
				require(present || p.termEquals("false", 10), "Require measured true/false same-generation model membership, not NA");
				increment(panel, member, present ? 3 : 2);
				final String full=p.parseString(9), id=CompareCobalaminBenchmarkGrades.modelId(full);
				Model model=byName.get(id);
				if(model==null){
					model=new Model(id, full, member, present, inputs.length);
					byName.put(id, model);
					models.add(model);
					counts[inputs.length][member][present ? 5 : 4]++;
					counts[inputs.length][3][present ? 5 : 4]++;
				}else{
					require(model.present==present && model.member==member && model.full.equals(full),
						"Original model has conflicting member, metadata or subset membership: "+id);
				}
				if(model.losses[panel]++==0){
					counts[panel][member][present ? 5 : 4]++;
					counts[panel][3][present ? 5 : 4]++;
				}
				if(present){presentOut.writer.print(buffer.clear().append(labels[panel]).tab().append(row).nl());}
			}
		}
	}

	/** Native fixed member slots, with ALL at slot three. */
	private int member(){
		for(int i=0; i<3; i++){if(p.termEquals(MEMBERS[i], 8)){return i;}}
		throw new IllegalArgumentException("Unexpected emitted member in baseline loss");
	}

	/** Adds a locus evidence row to member/panel and explicitly labeled union totals. */
	private void increment(final int panel, final int member, final int field){
		assert(field>=0 && field<4) : "Only locus fields are additive across panels; unique model counts are deduplicated separately";
		counts[panel][member][field]++;
		counts[panel][3][field]++;
		counts[inputs.length][member][field]++;
		counts[inputs.length][3][field]++;
	}

	private void write(){
		try(Output all=new Output(out+".models.tsv"); Output ids=new Output(out+".removed.ids.txt");
				Output main=new Output(out+".RF00174.ids.txt"); Output minor=new Output(out+".RF01482.ids.txt");
				Output third=new Output(out+".RF01689.ids.txt")){
			final Output[] memberIds={main, minor, third};
			all.writer.println("panel\temitted_member\tmodel_id\tno_overlap_lost_references\talready_in_subset\tstatus\tfull_original_model_text");
			for(final Model model:models){
				for(int i=0; i<inputs.length; i++){
					if(model.losses[i]>0){
						all.writer.print(buffer.clear().append(labels[i]).tab().append(MEMBERS[model.member]).tab().append(model.id)
							.tab().append(model.losses[i]).tab().append(model.present)
							.tab().append(model.present ? "ALREADY_PRESENT_NOT_AN_ADD_BACK" : "REMOVED_ORIGINAL_WINNER")
							.tab().append(model.full).nl());
					}
				}
				if(!model.present){
					ids.writer.println(model.id);
					memberIds[model.member].writer.println(model.id);
				}
			}
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("panel\temitted_member\tgraded_lost_reference_rows\tno_overlap_lost_reference_rows\tremoved_winner_rows\talready_present_winner_rows\tdistinct_removed_winners\tdistinct_already_present_winners\tqualification");
			for(int i=0; i<counts.length; i++){
				for(int member=0; member<4; member++){
					final long[] values=counts[i][member];
					assert(values[1]==values[2]+values[3]) : "Every absent-call loss is assigned exactly one measured subset-membership class";
					buffer.clear().append(i==inputs.length ? "UNION" : labels[i]).tab().append(MEMBERS[member]);
					for(final long value:values){buffer.tab().append(value);}
					summary.writer.print(buffer.append("\tREPORT_ONLY_UNION_LOCUS_COUNTS_ARE_EVIDENCE_ROWS\n"));
				}
			}
		}
	}

	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/** One retained identity with one primitive evidence count per source panel. */
	private static final class Model{
		Model(final String id_, final String full_, final int member_, final boolean present_, final int panels){
			id=id_;
			full=full_;
			member=member_;
			present=present_;
			losses=new long[panels];
		}
		final String id, full;
		final int member;
		final boolean present;
		final long[] losses;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String[] inputs, labels;
	private final String out;
	private final long[][][] counts;
	private final ObjectMap<String, Model> byName=new ObjectMap<String, Model>(String.class, Model.class);
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689", "ALL"};
	private static final String HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model";
}
