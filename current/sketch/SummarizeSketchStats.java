package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Colors;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import tax.TaxNode;
import tax.TaxTree;

/**
 * @author Brian Bushnell
 * @date June 28, 2017
 *
 */
public class SummarizeSketchStats {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new SummarizeSketchStats instance
		SummarizeSketchStats x=new SummarizeSketchStats(args);
		
		///And run it
		x.summarize();
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public SummarizeSketchStats(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		ArrayList<String> names=new ArrayList<String>();
		String taxTreeFile=null;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("printtotal") || a.equals("pt")){
				printTotal=Parse.parseBoolean(b);
			}else if(a.equals("ignoresametaxa")){
				ignoreSameTaxa=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamebarcode") || a.equals("ignoresameindex")){
				ignoreSameBarcode=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamelocation") || a.equals("ignoresameloc")){
				ignoreSameLocation=Parse.parseBoolean(b);
			}else if(a.equals("usetotal") || a.equals("totaldenominator") || a.equals("totald") || a.equals("td")){
				totalDenominator=Parse.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl") || a.equals("minlevel")){
				taxLevel=TaxTree.parseLevel(b);
				if(taxLevel>=0){
					taxLevel=TaxTree.levelToExtended(taxLevel);
				}
			}else if(a.equalsIgnoreCase("unique") || a.equalsIgnoreCase("uniquehits")){
				uniqueHitsForSecond=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("header") || a.equalsIgnoreCase("printheader")){
				printHeader=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){names.add(x2);}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		{//Process parser fields
			out=(parser.out1==null ? "stdout" : parser.out1);
			if(parser.in1!=null){
				String[] x=(new File(parser.in1).exists() ? new String[] {parser.in1} : parser.in1.split(","));
				for(String x2 : x){names.add(x2);}
			}
		}

		in=new ArrayList<String>();
		for(String s : names){
			Tools.getFileOrFiles(s, in, false, false, false, true);
		}
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
	}
	
	void setTaxtree(String taxTreeFile){
		if(taxTreeFile==null){
			return;
		}
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, false, false);
	}
	
	public void summarize(){
		ArrayList<SketchResultsSummary> list=new ArrayList<SketchResultsSummary>();
		for(String fname : in){
			ArrayList<SketchResultsSummary> ssl=summarize(fname);
			list.addAll(ssl);
		}
		
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		if(printHeader){tsw.print(header());}
//		if(printTotal){
//			tsw.println(total.toString());
//		}
		for(SketchResultsSummary ss : list){
			tsw.print(ss.toString());
		}
		tsw.poisonAndWait();
	}
	
//	Query: Troseus_1X_k55.fa	Seqs: 121 	Bases: 2410606	gSize: 2368581	SketchLen: 8923
//	WKID	KID	ANI	Complt	Contam	Matches	Unique	noHit	TaxID	gSize	gSeqs	taxName
//	99.89%	50.73%	100.00%	50.77%	0.02%	5683	5683	5	0	4719674	1	.	Troseus
	
	private ArrayList<SketchResultsSummary> summarize(String fname){
		TextFile tf=new TextFile(fname);
		ArrayList<SketchResultsSummary> list=new ArrayList<SketchResultsSummary>();
		SketchResultsSummary current=null;
		
		final String format="WKID	KID	ANI	Complt	Contam	Matches	Unique	noHit	TaxID	gSize	gSeqs	taxName";
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("Query:")){
				if(current!=null){list.add(current);}
				current=new SketchResultsSummary(line);
			}else if(line.startsWith("WKID")){
				assert(line.equals(format)) :
					"Format should be:\n"+format;
			}else if(line.length()>0){
				assert(current!=null) : "No Query Header for line "+line;
				current.add(line);
			}
		}
		if(current!=null){list.add(current);}
		tf.close();
		return list;
	}
	
	public static String header(){
		StringBuilder sb=new StringBuilder();
		
		sb.append("#query");

		sb.append('\t').append("seqs");
		sb.append('\t').append("bases");
		sb.append('\t').append("gSize");
		sb.append('\t').append("sketchLen");
		
		sb.append('\t').append("primaryHits");
		sb.append('\t').append("primaryUnique");
		sb.append('\t').append("primaryNoHit");

		sb.append('\t').append("WKID");
		sb.append('\t').append("KID");
		sb.append('\t').append("ANI");
		sb.append('\t').append("Complt");
		sb.append('\t').append("Contam");
		sb.append('\t').append("TaxID");
		sb.append('\t').append("TaxName");
		sb.append('\t').append("topContamID");
		sb.append('\t').append("topContamName");
		
		sb.append('\n');
		
		return sb.toString();
	}
	
	private class SketchResultsSummary {
		
		SketchResultsSummary(String line){
			parseHeader(line);
		}

		void parseHeader(String line){
			String[] split=line.split("\t");
			for(String s : split){
				String[] split2=s.trim().split(": ");
				assert(split2.length==2) : "\n"+line+"\n"+s+"\n"+Arrays.toString(split2)+"\n";
				String a=split2[0], b=split2[1];
//				outstream.println(a+", "+b);
				if(a.equals("Query")){
					query=b;
				}else if(a.equals("Seqs")){
					seqs=Integer.parseInt(b);
				}else if(a.equals("Bases")){
					bases=Long.parseLong(b);
				}else if(a.equals("gSize")){
					gSize=Long.parseLong(b);
				}else if(a.equals("SketchLen")){
					sketchLen=Integer.parseInt(b);
				}else if(a.equals("TaxID")){
					taxID=Integer.parseInt(b);
				}else if(a.equals("IMG")){
					img=Long.parseLong(b);
				}else if(a.equals("File")){
					sketchLen=Integer.parseInt(b);
				}
			}
		}
		
		public void add(String line) {
			SketchResultsLine srl=new SketchResultsLine(line);
			list.add(srl);
		}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append(query);

			sb.append('\t').append(seqs);
			sb.append('\t').append(bases);
			sb.append('\t').append(gSize);
			sb.append('\t').append(sketchLen);
			
			int primaryHits=0;
			int primaryUnique=0;
			int primaryNoHit=0;

			float WKID=0;
			float KID=0;
			float ANI=0;
			float Complt=0;
			float Contam=0;
			int TaxID=0;
			String TaxName=".";
			int topContamID=0;
			String topContamName=".";
			
			SketchResultsLine first=list.size()>0 ? list.get(0) : null;
			SketchResultsLine second=list.size()>1 ? list.get(1) : null;
			for(int i=2; tree!=null && i<list.size() && failsLevelFilter(first.taxID, second.taxID); i++){
				second=list.get(i);
			}
			if(second!=null && failsLevelFilter(first.taxID, second.taxID)){second=list.get(1);}
			
			if(second!=null && uniqueHitsForSecond){
				for(int i=1; i<list.size(); i++){
					
					SketchResultsLine line=list.get(i);
					if(!failsLevelFilter(first.taxID, line.taxID) && line.unique>second.unique && line.unique>=minUniqueHits){
						second=line;
					}
				}
			}
			
			if(first!=null){
				primaryHits=first.matches;
				primaryUnique=first.unique;
				primaryNoHit=first.noHit;

				WKID=first.wkid;
				KID=first.kid;
				ANI=first.ani;
				Complt=first.complt;
				Contam=first.contam;
				TaxID=first.taxID;
				TaxName=first.name;
			}
			if(second!=null){
				topContamID=second.taxID;
				topContamName=second.name;
			}
			
			sb.append('\t').append(primaryHits);
			sb.append('\t').append(primaryUnique);
			sb.append('\t').append(primaryNoHit);

			sb.append('\t').append(Tools.format("%.2f", WKID));
			sb.append('\t').append(Tools.format("%.2f", KID));
			sb.append('\t').append(Tools.format("%.2f", ANI));
			sb.append('\t').append(Tools.format("%.2f", Complt));
			sb.append('\t').append(Tools.format("%.2f", Contam));
			sb.append('\t').append(TaxID);
			sb.append('\t').append(TaxName);
			sb.append('\t').append(topContamID);
			sb.append('\t').append(topContamName);
			
			sb.append('\n');
			
			return sb.toString();
		}
		
		private boolean failsLevelFilter(int a, int b) {
			if(a<1 || b<1 || tree==null){return false;}
			int c=tree.commonAncestor(a, b);
			TaxNode tn=tree.getNode(c);
			while(!tn.cellularOrganisms() && tn.levelExtended==TaxTree.NO_RANK_E){tn=tree.getNode(tn.pid);}
			
			return tn.levelExtended<=taxLevel;
		}

		String query;
		String fname;
		int seqs;
		long bases;
		long gSize;
		int sketchLen;
		int taxID;
		long img;
		
		ArrayList<SketchResultsLine> list=new ArrayList<SketchResultsLine>();
		
	}
	
	private class SketchResultsLine{
		
		SketchResultsLine(String line){
			//Handle colors
			if(line.startsWith(Colors.esc)){
				int first=line.indexOf('m');
				int last=line.lastIndexOf(Colors.esc);
				line=line.substring(first+1, last);
			}
			String[] split=line.replaceAll("%", "").split("\t");
			wkid=Float.parseFloat(split[0]);
			kid=Float.parseFloat(split[1]);
			ani=Float.parseFloat(split[2]);
			complt=Float.parseFloat(split[3]);
			contam=Float.parseFloat(split[4]);
			
			matches=Integer.parseInt(split[5]);
			unique=Integer.parseInt(split[6]);
			noHit=Integer.parseInt(split[7]);
			taxID=Integer.parseInt(split[8]);
			gSize=Integer.parseInt(split[9]);
			gSeqs=Integer.parseInt(split[10]);
			
			name=split[11];
			if(name.equals(".") && split.length>11){
				name=split[12];
			}
		}
		
		float wkid;
		float kid;
		float ani;
		float complt;
		float contam;
		int matches;
		int unique;
		int noHit;
		int taxID;
		int gSize;
		int gSeqs;
		String name;
	}
	
	final ArrayList<String> in;
	final String out;
	
	TaxTree tree=null;
	int taxLevel=TaxTree.GENUS_E;
	boolean uniqueHitsForSecond=false;
	int minUniqueHits=3;
	boolean printHeader=true;
	
	/** Legacy code from SealStats */
	boolean ignoreSameTaxa=false;
	boolean ignoreSameBarcode=false;
	boolean ignoreSameLocation=false;
	boolean totalDenominator=false;
	boolean printTotal=true;
	
	PrintStream outstream=System.err;
	
}
