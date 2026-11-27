package prok;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;

/** Immutable static endpoint routes on one complete native consensus axis.
 * The table is an already selected artifact; loading never chooses a route.
 * @author Brian, Ganyu
 */
public final class NcrnaModelRoutes{

	/*--------------------------------------------------------------*/
	/*----------------           Loading            ----------------*/
	/*--------------------------------------------------------------*/

	public static NcrnaModelRoutes load(final String path, final String member, final String[] names){
		require(path!=null && ("RF00174".equals(member) || "RF01482".equals(member)) && names!=null && names.length>0,
			"Static hybrid routes require one declared refined member and its complete ordered consensus axis");
		final boolean[] tip=new boolean[names.length];
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder encoded=new ByteBuilder();
		int count=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, HEADER, path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==11 && (p.termEquals("RF00174", 0) || p.termEquals("RF01482", 0)),
					"Frozen route rows must retain their complete schema and one of the two refined member identities");
				if(!p.termEquals(member, 0)){continue;}
				require(count<names.length && p.parseInt(1)==count
					&& p.termEquals(BuildNcrnaReference.field(encoded.clear(), NcrnaSupport.recordName(names[count])).toString(), 2)
					&& (p.termEquals("TIP", 3) || p.termEquals("CM", 3)),
					"Every native model must have exactly one frozen route in its original order: member="+member+", ordinal="+count);
				tip[count++]=p.termEquals("TIP", 3);
			}
		}
		require(count==names.length, "Frozen routes must include sparse and unobserved models, not just observed winners");
		return new NcrnaModelRoutes(tip);
	}
	private NcrnaModelRoutes(final boolean[] tip_){tip=tip_;}

	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/

	public boolean tip(final int model){
		require(model>=0 && model<tip.length, "Call-time and teaching routes must name a real model on the frozen axis");
		return tip[model];
	}
	public int size(){return tip.length;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean[] tip;
	public static final String HEADER="member\tmodel_ordinal0\tmodel\troute\treason\tmatched\tresolved_groups\tcandidates\tselected\ttip_cpu_ns\tcm_cpu_ns";
}
