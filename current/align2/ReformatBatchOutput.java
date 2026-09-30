package align2;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextFile;
import shared.Tools;

/**
 * Reformats legacy GradeSamFile batch reports into a tab-separated table.
 * Expects labelled tab-delimited metrics and experiment
 * names with underscore-separated program, mutation, and read-count fields.
 * This is not a general parser for current mapper logs or arbitrary filenames.
 * Complete blocks are emitted as 17-column rows; absent optional strict/loose
 * correctness metrics remain empty fields. Incomplete or duplicate fields fail
 * before their row is emitted. Input is materialized in memory; output goes to stdout.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class ReformatBatchOutput{

	public static void main(String[] args){
		assert(args.length>0) : "A legacy batch-report input path is required";
		TextFile tf=new TextFile(args[0], false);
		final String[] lines;
		try{lines=tf.toStringLines();}
		finally{tf.close();}
		ArrayList<String> list=new ArrayList<String>();

		int mode=0;

		System.out.println(header());

		for(String s : lines){
			if(mode==0 && s.startsWith("Mapping Statistics for ")){
				throw new IllegalArgumentException("Missing Elapsed line before batch-report name: "+s);
			}
			if(s.startsWith("Elapsed:")){
				if(!list.isEmpty()){
					process(list); //Reject an incomplete previous block before starting another.
					list.clear();
					mode=0;
				}
				mode++;
			}

			if(mode>0){
				list.add(s);
				if(s.startsWith("false negative:")){
					process(list);
					list.clear();
					mode=0;
				}
			}
		}
		if(!list.isEmpty()){process(list);}
	}

	public static String header(){
		return("program\tfile\tvartype\tcount\treads\tprimary\tsecondary\ttime\tmapped\tretained\tdiscarded\tambiguous\ttruePositive\t" +
				"falsePositive\ttruePositiveL\tfalsePositiveL\tfalseNegative");
	}

	/**
	 * Extracts read count from a legacy experiment name after removal of .sam:.
	 * Accepts a final rCOUNT or rCOUNTxLENGTH field, or a COUNTxLENGTHbp field.
	 * For example, bwa_1S_0I_r400000x100 has 400000 reads.
	 *
	 * @param name Filename to parse
	 * @return Number of reads, or 0 if not found
	 */
	public static int getReads(String name){
		String[] split=name.split("_");
		String r=(split[split.length-1]);
		if(r.charAt(0)=='r' && Tools.isDigit(r.charAt(r.length()-1))){
			assert(r.charAt(0)=='r') : Arrays.toString(split)+", "+name;
			r=r.substring(1);
			if(r.contains("x")){
				r=r.substring(0, r.indexOf('x'));
			}
			return Integer.parseInt(r);
		}else{
			for(String s : split){
				if(s.endsWith("bp") && s.contains("x") && Tools.isDigit(s.charAt(0))){
					r=s.substring(0, s.indexOf('x'));
					return Integer.parseInt(r);
				}
			}
		}
		return 0;
	}

	/**
	 * Extracts variant type from filename based on BBTools naming convention.
	 * Identifies the variant type character (S, I, D, U, N) from filename components
	 * that start with non-zero digits.
	 *
	 * @param name Filename to parse
	 * @return Variant type character, or '?' if not found
	 */
	public static char getVarType(String name){
		//TODO: Probable bug - empty fields throw, and arbitrary nonzero numeric fields are treated as variants.
		//getCount shares this historical grammar; tightening it needs real legacy naming fixtures.
		String[] split=name.split("_");
		for(String s : split){
			char c=s.charAt(0);
			if(Tools.isDigit(c) && c!='0' && !s.endsWith("bp")){
				return s.charAt(s.length()-1);
			}
		}
		return '?';
	}

	/**
	 * Extracts variant count from filename based on BBTools naming convention.
	 * Parses numerical prefix from filename components to determine variant count.
	 * @param name Filename to parse
	 * @return Variant count, or 0 if not found
	 */
	public static int getCount(String name){
		String[] split=name.split("_");
		for(String s : split){
			char c=s.charAt(0);
			if(Tools.isDigit(c) && c!='0' && !s.endsWith("bp")){
				String r=s.substring(0, s.length()-1);
				return Integer.parseInt(r);
			}
		}
		return 0;
	}

	/**
	 * Extracts program name from filename.
	 * Returns substring before the first underscore character.
	 * @param name Filename to parse
	 * @return Program name prefix
	 */
	public static String getProgram(String name){
		return name.substring(0, name.indexOf('_'));
	}

	/**
	 * Processes a block of mapping statistics and outputs formatted results.
	 * Parses elapsed time, alignment counts, mapping percentages, and correctness
	 * metrics from the input lines and outputs a tab-delimited summary.
	 * @param list List of lines containing mapping statistics for one alignment run
	 */
	public static void process(ArrayList<String> list){

		String name=null;
		String time=null;
		final String[] metrics=new String[9];
		int primary=0, secondary=0, correctnessFields=0;
		boolean havePrimary=false, haveSecondary=false;
		boolean haveStrictHeader=false, haveLooseHeader=false;

		for(String s : list){
			String[] split=s.split("\t", -1);
			if(split.length>2){throw new IllegalArgumentException("Too many tab-delimited fields in batch-report line: "+s);}
			String a=split[0].trim();
			String b=split.length>1 ? split[1].trim() : null;
			if(a.equals("Elapsed:")){
				if(time!=null){throw new IllegalArgumentException("Duplicate Elapsed field in batch report");}
				time=requiredValue(b, s);
			}else if(a.startsWith("lines:") || a.equals("Mapping:")){
				//do nothing
			}else if(a.startsWith("Strict correctness")){
				if(haveStrictHeader || correctnessFields!=0){throw new IllegalArgumentException("Misplaced strict correctness heading: "+s);}
				haveStrictHeader=true;
			}else if(a.startsWith("Loose correctness")){
				if(haveLooseHeader || correctnessFields!=2){throw new IllegalArgumentException("Misplaced loose correctness heading: "+s);}
				haveLooseHeader=true;
			}else if(a.startsWith("Mapping Statistics for ")){
				if(b!=null){throw new IllegalArgumentException("Tab in batch-report experiment name: "+s);}
				if(name!=null){throw new IllegalArgumentException("Duplicate experiment name in batch report: "+s);}
				name=a.replace("Mapping Statistics for ", "").replace(".sam:", "");
			}else if(a.startsWith("primary alignments:")){
				if(havePrimary){throw new IllegalArgumentException("Duplicate primary alignment count: "+s);}
				b=requiredValue(b, s);
				b=b.replace(" found of ", "_");
				b=b.replace(" expected", "");
				String[] split2=b.split("_");
				if(split2.length!=2){throw new IllegalArgumentException("Expected primary count and expected-read count: "+s);}
				primary=Integer.parseInt(split2[0]);
				final int expected=Integer.parseInt(split2[1]);
				if(primary<0 || expected<0){throw new IllegalArgumentException("Negative alignment/read count: "+s);}
				havePrimary=true;
			}else if(a.startsWith("secondary alignments:")){
				if(haveSecondary){throw new IllegalArgumentException("Duplicate secondary alignment count: "+s);}
				b=requiredValue(b, s);
				b=b.replace(" found", "");
				secondary=Integer.parseInt(b);
				if(secondary<0){throw new IllegalArgumentException("Negative secondary alignment count: "+s);}
				haveSecondary=true;
			}else{
				int field=-1;
				if(a.equals("mapped:")){field=0;}
				else if(a.equals("retained:")){field=1;}
				else if(a.equals("discarded:")){field=2;}
				else if(a.equals("ambiguous:")){field=3;}
				else if(a.equals("false negative:")){field=8;}
				else if(a.equals("true positive:") || a.equals("false positive:")){
					final String expectedLabel=(correctnessFields%2==0 ? "true positive:" : "false positive:");
					if(correctnessFields>=4 || !a.equals(expectedLabel)){
						throw new IllegalArgumentException("Expected strict then loose true/false-positive pairs: "+s);
					}
					field=4+correctnessFields++;
				}
				if(field>=0){
					if(metrics[field]!=null){throw new IllegalArgumentException("Duplicate batch-report metric: "+s);}
					metrics[field]=requiredValue(requiredValue(b, s).replace("%", "").trim(), s);
				}else if(b!=null){throw new IllegalArgumentException("Unrecognized tab-delimited batch-report field: "+s);}
			}

		}
		if(name==null || time==null || !havePrimary || !haveSecondary){
			throw new IllegalArgumentException("Incomplete batch report: expected name, elapsed time, primary and secondary counts; name="+name);
		}
		for(int i : new int[]{0, 1, 2, 3, 8}){
			if(metrics[i]==null){throw new IllegalArgumentException("Incomplete batch report for "+name+": missing "+METRIC_NAMES[i]);}
		}
		//GradeSamFile emits all four truth metrics only when parsecustom is enabled.
		if(correctnessFields!=4 && (correctnessFields!=0 || haveStrictHeader || haveLooseHeader)){
			throw new IllegalArgumentException("Incomplete strict/loose correctness section for "+name);
		}

		String prg=null;
		char type='S';
		int reads=1;
		int vars=0;

		if(name!=null){
			try{
				prg=getProgram(name);
				type=getVarType(name);
				reads=getReads(name);
				vars=getCount(name);
			}catch(RuntimeException e){
				throw new IllegalArgumentException("Cannot parse legacy batch-report experiment name: "+name, e);
			}
		}

		final StringBuilder sb=new StringBuilder();
		for(String metric : metrics){sb.append('\t'); if(metric!=null){sb.append(metric);}}
		System.out.println(prg+"\t"+name+"\t"+type+"\t"+vars+"\t"+reads+"\t"+primary+"\t"+secondary+"\t"+time+sb);

	}

	/** A labelled field cannot silently become an empty table cell. */
	private static String requiredValue(final String value, final String line){
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing batch-report value: "+line);}
		return value;
	}

	private static final String[] METRIC_NAMES={"mapped", "retained", "discarded", "ambiguous",
			"truePositive", "falsePositive", "truePositiveL", "falsePositiveL", "falseNegative"};

}
