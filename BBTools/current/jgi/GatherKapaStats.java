    package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import json.JsonObject;
import json.JsonParser;
import server.ServerTools;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class GatherKapaStats {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		GatherKapaStats x=new GatherKapaStats(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public GatherKapaStats(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				ByteFile1.verbose=verbose;
//				ByteFile2.verbose=verbose;
//				ReadWrite.verbose=verbose;
			}else if(a.equals("raw") || a.equals("printraw")){
				printRaw=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ArrayList<Plate> plates=loadPlates(bf);
		errorState|=bf.close();
		
		analyzePlates(plates);
		
		ByteStreamWriter bsw=makeBSW(ffout1);
		if(bsw!=null){
			if(printRaw){
				printRawResults(bsw);
			}else{
				printResults(bsw);
			}
			errorState|=bsw.poisonAndWait();
		}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Lines Out:         \t"+linesOut);
//		outstream.println("Valid Lines:       \t"+linesOut);
//		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<Plate> loadPlates(ByteFile bf){
		
		ArrayList<Plate> plates=new ArrayList<Plate>();
		byte[] line=bf.nextLine();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				
				if(valid){
					String[] split=new String(line).split("\t");
					assert(split.length>=1) : Arrays.toString(split);
					String name=split[0];
					String lot=(split.length>1 ? split[1] : null);
					Plate plate=new Plate(name, lot);
					plate.fillFromWeb();
					plates.add(plate);
					plateMap.put(name, plate);
				}
			}
			line=bf.nextLine();
		}
		return plates;
	}
	
	private void analyzePlates(ArrayList<Plate> plates){
		for(Plate p : plates){
			for(Well w : p.wells){
				long kapaReads=w.correctKapaReads+w.incorrectKapaReads;
				if(kapaReads>0){
					final double mult=1000000.0/kapaReads;
					TagData td=tagMap.get(w.correctKapaTag);
					if(td==null){
						td=new TagData(w.correctKapaTag, w.name);
						tagMap.put(w.correctKapaTag, td);
					}
					td.timesSeen++;
					for(Entry<String, KapaEntry> e : w.kapaMap.entrySet()){
						KapaEntry ke=e.getValue();
						double ppmk=ke.reads*mult;
						td.add(ke.tagName, ppmk, p.name);
					}
				}
			}
		}
	}
	
	private void printResults(ByteStreamWriter bsw){
		ArrayList<TagData> list=new ArrayList<TagData>();
		list.addAll(tagMap.values());
		Collections.sort(list);
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Tag\tOther\tMin\t25%\t50%\t75%\tMax\tAvg\tStdev\tObserved\tTotal\tFraction\n");
		for(TagData td : list){
			ArrayList<String> keys=new ArrayList<String>();
			keys.addAll(td.ppmMap.keySet());
			Collections.sort(keys);
			final int len=td.timesSeen;
			for(String key : keys){
				double[] ppmk=td.getPpmArray(key, true);
//				if(ppmk.length<seen){
//					ppmk=Arrays.copyOf(ppmk, newLength)
//				}
				assert(len==ppmk.length);
				int count=0;
				for(double d : ppmk){
					if(d>0){count++;}
				}
				double min=ppmk[0];
				double max=ppmk[len-1];
				double p25=ppmk[(int)Math.round((len-1)*.25)];
				double p50=ppmk[(int)Math.round((len-1)*.50)];
				double p75=ppmk[(int)Math.round((len-1)*.75)];
				double avg=shared.Vector.sum(ppmk)/len;
				double stdev=Tools.standardDeviation(ppmk);
				bb.append(td.name).append('\t');
				bb.append(key).append('\t');
				bb.append(min, 2).append('\t');
				bb.append(p25, 2).append('\t');
				bb.append(p50, 2).append('\t');
				bb.append(p75, 2).append('\t');
				bb.append(max, 2).append('\t');
				bb.append(avg, 2).append('\t');
				bb.append(stdev, 2).append('\t');
				bb.append(count).append('\t');
				bb.append(len).append('\t');
				bb.append(count/(double)len, 4).append('\n');
				bsw.print(bb);
				linesOut++;
				bytesOut+=bb.length;
				bb.clear();
//				bsw.println(Arrays.toString(ppmk));
			}
		}
		if(!bb.isEmpty()){
			linesOut++;
			bytesOut+=bb.length;
			bsw.print(bb);
		}
	}
	
	private void printRawResults0(ByteStreamWriter bsw){
		ArrayList<TagData> list=new ArrayList<TagData>();
		list.addAll(tagMap.values());
		Collections.sort(list);
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Tag\tOther\tTotal\tPPM,...\n");
		for(TagData td : list){
			ArrayList<String> keys=new ArrayList<String>();
			keys.addAll(td.ppmMap.keySet());
			Collections.sort(keys);
			final int len=td.timesSeen;
			for(String key : keys){
				double[] ppmk=td.getPpmArray(key, true);
//				if(ppmk.length<seen){
//					ppmk=Arrays.copyOf(ppmk, newLength)
//				}
				assert(len==ppmk.length);
				bb.append(td.name).append('\t');
				bb.append(key).append('\t');
				bb.append(len).append('\t');
				String comma="";
				for(double d : ppmk){
					bb.append(comma);
					bb.append(d, 2);
					comma=",";
				}
				bb.append('\n');
				
				bsw.print(bb);
				linesOut++;
				bytesOut+=bb.length;
				bb.clear();
//				bsw.println(Arrays.toString(ppmk));
			}
		}
		if(!bb.isEmpty()){
			linesOut++;
			bytesOut+=bb.length;
			bsw.print(bb);
		}
	}
	
	private void printRawResults(ByteStreamWriter bsw){
		ArrayList<TagData> list=new ArrayList<TagData>();
		list.addAll(tagMap.values());
		Collections.sort(list);
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Plate\tSinkWell\tSinkCorrectTag\tSinkReads\tSinkCorrectKapaReads\tSinkTotalKapaReads\t"
				+ "SourceWell\tMeasuredTag\tSourceReads\tSourceCorrectKapaReads\tSourceKapaReadsInSink\t"
				+ "KPPM (SourceKapa/SinkKapa)\tGReads (InferredContamGenomicReads)\t"
				+ "GPPM (InferredContamGenomicPPM)\n");
		for(TagData td : list){
			ArrayList<String> keys=new ArrayList<String>();
			keys.addAll(td.ppmMap.keySet());
			Collections.sort(keys);
			final int len=td.timesSeen;
			for(String key : keys){
				double[] ppmk=td.getPpmArray(key, false);
				String[] plateNames=td.getPlateNameArray(key, false);
//				if(ppmk.length<seen){
//					ppmk=Arrays.copyOf(ppmk, newLength)
//				}
				assert(len==ppmk.length);
				for(int i=0; i<ppmk.length; i++){
					double d=ppmk[i];
					String plateName=plateNames[i];
					if(plateName!=null){
						Plate plate=plateMap.get(plateName);
						Well sink=plate.tagToCorrectWellMap.get(td.name);
						Well source=plate.tagToCorrectWellMap.get(key);
						if(source==null){source=dummy;}
						KapaEntry keSource=sink.kapaMap.get(key);
						long contamReads=keSource.reads;
						double greads=contamReads*(source.reads/(double)source.correctKapaReads);
						double gppm=1000000*greads/(double)sink.reads;
						
						bb.append(plateName).tab();
						bb.append(td.wellName).tab();
						bb.append(td.name).tab();
						bb.append(sink.reads).tab();
						bb.append(sink.correctKapaReads).tab();
						bb.append(sink.correctKapaReads+sink.incorrectKapaReads).tab();
						bb.append(source.name).tab();
						bb.append(key).tab();
						bb.append(source.reads).tab();
						bb.append(source.correctKapaReads).tab();
						bb.append(contamReads).tab();
						bb.append(d, 2).tab();
						bb.append(greads, 2).tab();
						bb.append(gppm, 2).nl();
					}
				}
				
				bsw.print(bb);
				linesOut++;
				bytesOut+=bb.length;
				bb.clear();
//				bsw.println(Arrays.toString(ppmk));
			}
		}
		if(!bb.isEmpty()){
			linesOut++;
			bytesOut+=bb.length;
			bsw.print(bb);
		}
	}
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	class Plate{
		
		public Plate(String name_, String lot_){
			name=name_;
			lot=lot_;
		}
		
		void fillFromWeb(){
			JsonObject data=grabData();
			int size=data.jmapSize();
			wells=new ArrayList<Well>(size);
			if(size<1){
				outstream.println("No Kapa for plate "+name);
				return;
			}
			for(Entry<String, JsonObject> e : data.jmap.entrySet()){
				String key=e.getKey();
				JsonObject jo=e.getValue();
				Well well=new Well(key, jo, name);
				wells.add(well);
				tagToCorrectWellMap.put(well.correctKapaTag, well);
				if(verbose && well.name.contentEquals("B1")){
					System.err.println(well);
				}
			}
		}
		
		JsonObject grabData(){
			String address=addressPrefix+name+addressSuffix;
//			System.err.println("Reading "+address);
			ByteBuilder message=ServerTools.readPage(address, true);
			assert(message!=null && message.length()>0) : "No data from address "+address;
			
			jp.set(message.toBytes());
			JsonObject jo=jp.parseJsonObject();
			assert(jo!=null && jo.jmapSize()==1 && jo.omapSize()==0) : jo.toString();
			
			JsonObject data=jo.getJson("data");
			assert(data!=null) : jo.toString();
//			assert(data.jmapSize()>0 /*&& data.smapSize()==0*/) : data.toString()+"\n\n"+data.smap+"\n\n"+data.jmapSize(); //These assertions are not important, just making sure I understand the API
			return data;
		}

		final String name;
		final String lot;
		
		ArrayList<Well> wells;
		LinkedHashMap<String, Well> tagToCorrectWellMap=new LinkedHashMap<String, Well>();
		
	}
	
	class Well{
		
//		"library_name":"CSNGW",
//        "asm_comments":"",
//        "instrument_type":"HiSeq-2500 1TB",
//        "asm_qc_status":"Pass",
//        "dt_created":"2018-05-08 18:05:36",
//        "alq_container_barcode":"27-353939",
//        "seq_unit_name":"12396.3.255039.ATGCCTG-ACAGGCA.fastq.gz",
//        "seq_proj_id":"1190410",
//        "raw_reads":8039600,
//        "seq_proj_name":"Ensifer meliloti 417 Resequencing",
//        "account_jgi_sci_prog":"Microbial",
//        "alq_initial_mass_ng":"n.a.",
//        "library_protocol":"Regular (DNA)",
//        "run_configuration":"2x101",
		
		public Well(String name_, JsonObject jo, String plate){
			name=name_;

			library=jo.getString("library_name");
			instrument=jo.getString("instrument_type");
			date=jo.getString("dt_created");
			alq_container_barcode=jo.getString("alq_container_barcode");
			seq_unit_name=jo.getString("seq_unit_name");
			seq_proj_id=jo.getString("seq_proj_id");
			seq_proj_name=jo.getString("seq_proj_name");
			Long temp=jo.getLong("raw_reads");
			reads=(temp==null ? 0 : temp.longValue());
			run_configuration=jo.getString("run_configuration");
			
			if(name.equalsIgnoreCase("X")){return;}
			JsonObject kapa=jo.getJson("kapa");
			if(kapa==null && outstream!=null){
				outstream.println("No Kapa for "+library+", plate "+plate);
			}else{
				loadKapa(kapa);
			}
		}
		
//		"hit":77972,
//        "name":"tag059",
//        "offppm":3.9802975272401615,
//        "offhit":32,
//        "converted_offtarget_reads_ppm":246.31785207079872,
//        "kapa_stats_file":"dna/00/50/16/18//29495471-kapa.stats",
//        "pct":0.9698492462311559
		
		void loadKapa(JsonObject kapa){
			correctKapaTag=kapa.getString("name");
			correctKapaReads=kapa.getLong("hit");
			incorrectKapaReads=kapa.getLong("offhit");
			Number n=kapa.getNumber("converted_offtarget_reads_ppm");
			Class<?> c=n.getClass();
			if(c==Double.class){
				converted_offtarget_reads_ppm=(Double)n;
			}else{
				converted_offtarget_reads_ppm=(Long)n;
			}
			Object[] offwells=kapa.getArray("offwells");
			
			kapaMap=new LinkedHashMap<String, KapaEntry>(3+offwells.length*2);
			kapaMap.put(correctKapaTag, new KapaEntry(name, correctKapaReads, correctKapaTag));
//			tagToReads.put(correctKapaTag, correctKapaReads);
			for(Object o : offwells){
				KapaEntry ke=new KapaEntry((Object[])o);
//				tagToReads.put(ke.tagName, ke.reads);
				kapaMap.put(ke.tagName, ke);
			}
		}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder();
			sb.append("name\t"+name).append('\n');
			sb.append("correctKapaTag\t"+correctKapaTag).append('\n');
			sb.append("reads\t"+reads).append('\n');
			sb.append("correctKapaReads\t"+correctKapaReads).append('\n');
			sb.append("incorrectKapaReads\t"+incorrectKapaReads).append('\n');
			for(KapaEntry e : kapaMap.values()){
				sb.append(e.toString()).append('\n');
			}
			return sb.toString();
		}
		
		final String name;
		
		String library;
		String instrument;
		String date;
		String alq_container_barcode;
		String seq_unit_name;
		String seq_proj_id;
		long reads;
		String seq_proj_name;
		String run_configuration;
		
		String correctKapaTag;
		long correctKapaReads;
		long incorrectKapaReads;
		double converted_offtarget_reads_ppm;

		LinkedHashMap<String, KapaEntry> kapaMap;
//		LinkedHashMap<String, Long> tagToReads=new LinkedHashMap<String, Long>();
		
//        "asm_comments":"",
//        "instrument_type":"HiSeq-2500 1TB",
//        "asm_qc_status":"Pass",
//        "dt_created":"2018-05-08 18:05:36",
//        "alq_container_barcode":"27-353939",
//        "seq_unit_name":"12396.3.255039.ATGCCTG-ACAGGCA.fastq.gz",
//        "seq_proj_id":"1190410",
//        "raw_reads":8039600,
//        "seq_proj_name":"Ensifer meliloti 417 Resequencing",
//        "account_jgi_sci_prog":"Microbial",
//        "alq_initial_mass_ng":"n.a.",
//        "library_protocol":"Regular (DNA)",
//        "run_configuration":"2x101",
		
	}
	
	class KapaEntry{
		
//		"offwells":[
//                    [
//                        "A1",
//                        0.00017413801681675706,
//                        14,
//                        "tag001"
//                    ],

		KapaEntry(Object[] array){
			assert(array.length==4) : Arrays.toString(array);
			wellName=(String)array[0];
			reads=(Long)array[2];
			tagName=(String)array[3];
		}
		
		KapaEntry(String wellName_, long reads_, String tagName_){
			wellName=wellName_;
			reads=reads_;
			tagName=tagName_;
		}
		
		@Override
		public String toString(){
			return wellName+"\t"+tagName+"\t"+reads;
		}
		
		String wellName;
		long reads;
		String tagName;
		
	}
	
	//Ugly because it was retrofitted to support unsorted arrays and plate names
	class TagData implements Comparable<TagData> {
		
		TagData(String name_, String wellName_){
			name=name_;
			wellName=wellName_;
		}
		
		void add(String tag, double ppmk, String plate){
			ArrayList<Double> list=ppmMap.get(tag);
			if(list==null){
				list=new ArrayList<Double>();
				ppmMap.put(tag, list);
			}
			list.add(ppmk);

			ArrayList<String> list2=plateNameMap.get(tag);
			if(list2==null){
				list2=new ArrayList<String>();
				plateNameMap.put(tag, list2);
			}
			list2.add(plate);
		}
		
		double[] getPpmArray(String key, boolean sort){
			ArrayList<Double> list=ppmMap.get(key);
			return toPpmArray(list, sort);
		}
		
		String[] getPlateNameArray(String key, boolean sort){
			ArrayList<String> list=plateNameMap.get(key);
			return toPlateArray(list, sort);
		}
		
		double[] toPpmArray(ArrayList<Double> list, boolean sort){
			if(list==null){return null;}
//			double[] array=new double[list.size()];
			double[] array=new double[timesSeen];
			for(int i=0; i<list.size(); i++){
				array[i]=list.get(i);
			}
			if(sort){Arrays.sort(array);}
			return array;
		}
		
		String[] toPlateArray(ArrayList<String> list, boolean sort){
			if(list==null){return null;}
//			double[] array=new double[list.size()];
			String[] array=new String[timesSeen];
			for(int i=0; i<list.size(); i++){
				array[i]=list.get(i);
			}
			if(sort){Arrays.sort(array);}
			return array;
		}

		@Override
		public int compareTo(TagData other) {
			return name.compareTo(other.name);
		}
		
		final String name;
		final String wellName;
		int timesSeen=0;

		LinkedHashMap<String, ArrayList<Double>> ppmMap=new LinkedHashMap<String, ArrayList<Double>>(203);
		LinkedHashMap<String, ArrayList<String>> plateNameMap=new LinkedHashMap<String, ArrayList<String>>(203);
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private String addressPrefix="https://rqc.jgi.doe.gov/api/plate_ui/page/";
	private String addressSuffix="/kapa spikein";//"/kapa spikein";
	
	private boolean printRaw=false;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	private final JsonParser jp=new JsonParser();

	private final LinkedHashMap<String, TagData> tagMap=new LinkedHashMap<String, TagData>(203);
	private final LinkedHashMap<String, Plate> plateMap=new LinkedHashMap<String, Plate>(203);
	
	final Well dummy=new Well("X", new JsonObject(), "X");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
