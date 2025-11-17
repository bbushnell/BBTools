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
	
	/**
	 * Main processing pipeline for Kapa contamination analysis.
	 * Loads plate data from input, analyzes contamination patterns,
	 * and writes results to output file with timing information.
	 * @param t Timer for tracking execution time
	 */
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
	
	/**
	 * Loads plate information from input file and fetches detailed data from web API.
	 * Parses tab-separated input with plate names and lot numbers,
	 * then fetches full Kapa data from JGI web service.
	 *
	 * @param bf ByteFile containing plate names and lot numbers
	 * @return ArrayList of populated Plate objects with well data
	 */
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
	
	/**
	 * Analyzes contamination patterns across all plates and wells.
	 * Calculates parts-per-million contamination rates between different
	 * Kapa tags and builds statistical data for reporting.
	 * @param plates List of plates with well data to analyze
	 */
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
	
	/**
	 * Prints statistical summary of contamination results.
	 * Outputs quartile statistics, averages, and standard deviations
	 * for each tag-to-tag contamination relationship.
	 * @param bsw ByteStreamWriter for output
	 */
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
	
	/**
	 * Prints raw contamination data in comma-separated format.
	 * Alternative raw output format showing all individual measurements
	 * for each tag pair.
	 * @param bsw ByteStreamWriter for output
	 */
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
	
	/**
	 * Prints detailed raw results with full well information.
	 * Includes plate names, well identifiers, read counts, and calculated
	 * contamination metrics in tabular format for detailed analysis.
	 * @param bsw ByteStreamWriter for output
	 */
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
	
	/**
	 * Creates and starts a ByteStreamWriter for the specified file format.
	 * @param ff FileFormat specification for output
	 * @return Started ByteStreamWriter or null if ff is null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Represents a sequencing plate with associated wells and Kapa spike-in data.
	 * Fetches detailed well information from JGI web API and maintains
	 * mappings between Kapa tags and wells for contamination analysis.
	 */
	class Plate{
		
		/**
		 * Constructs a Plate with name and lot number.
		 * @param name_ Plate identifier
		 * @param lot_ Lot number for the plate
		 */
		public Plate(String name_, String lot_){
			name=name_;
			lot=lot_;
		}
		
		/**
		 * Fetches complete well data from JGI web API and populates plate information.
		 * Creates Well objects for each position and builds tag-to-well mapping
		 * for contamination analysis.
		 */
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
		
		/**
		 * Retrieves JSON data from JGI web API for this plate.
		 * Constructs API URL, fetches data, and parses JSON response
		 * containing well and Kapa information.
		 * @return JsonObject containing plate data from web API
		 */
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

		/** Plate identifier name */
		final String name;
		/** Plate lot number */
		final String lot;
		
		/** List of wells on this plate */
		ArrayList<Well> wells;
		/** Maps Kapa tags to their corresponding correct wells */
		LinkedHashMap<String, Well> tagToCorrectWellMap=new LinkedHashMap<String, Well>();
		
	}
	
	/**
	 * Represents a single well on a sequencing plate with library and Kapa data.
	 * Contains sequencing metadata, read counts, and Kapa contamination information
	 * including cross-contamination measurements from other wells.
	 */
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
		
		/**
		 * Constructs a Well from JSON data retrieved from JGI API.
		 * Parses library metadata, sequencing information, and Kapa contamination data.
		 *
		 * @param name_ Well position identifier (e.g., "A1", "B2")
		 * @param jo JsonObject containing well data from API
		 * @param plate Plate name for error reporting
		 */
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
		
		/**
		 * Loads Kapa spike-in contamination data from JSON object.
		 * Parses correct and incorrect Kapa reads, contamination rates,
		 * and cross-contamination data from other wells.
		 * @param kapa JsonObject containing Kapa contamination statistics
		 */
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
		
		/** Well position identifier */
		final String name;
		
		/** Library name for this well */
		String library;
		/** Sequencing instrument type used */
		String instrument;
		/** Date when sequencing was performed */
		String date;
		/** Container barcode for the well */
		String alq_container_barcode;
		/** Sequencing unit name/filename */
		String seq_unit_name;
		/** Sequencing project identifier */
		String seq_proj_id;
		/** Total number of raw reads for this well */
		long reads;
		/** Sequencing project name */
		String seq_proj_name;
		/** Sequencing run configuration (e.g., "2x101") */
		String run_configuration;
		
		/** The correct Kapa tag assigned to this well */
		String correctKapaTag;
		/** Number of reads with correct Kapa tag */
		long correctKapaReads;
		/** Number of reads with incorrect Kapa tags */
		long incorrectKapaReads;
		/** Parts per million of converted off-target reads */
		double converted_offtarget_reads_ppm;

		/** Maps Kapa tag names to contamination entries for this well */
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
	
	/**
	 * Represents contamination data for a specific Kapa tag in a well.
	 * Contains the source well name, number of contaminating reads,
	 * and the tag identifier for cross-contamination analysis.
	 */
	class KapaEntry{
		
//		"offwells":[
//                    [
//                        "A1",
//                        0.00017413801681675706,
//                        14,
//                        "tag001"
//                    ],

		/**
		 * Constructs KapaEntry from JSON array data.
		 * Parses well name, read count, and tag name from API response format.
		 * @param array Object array containing [wellName, ?, readCount, tagName]
		 */
		KapaEntry(Object[] array){
			assert(array.length==4) : Arrays.toString(array);
			wellName=(String)array[0];
			reads=(Long)array[2];
			tagName=(String)array[3];
		}
		
		/**
		 * Constructs KapaEntry with explicit parameters.
		 * @param wellName_ Source well name for contamination
		 * @param reads_ Number of contaminating reads
		 * @param tagName_ Kapa tag identifier
		 */
		KapaEntry(String wellName_, long reads_, String tagName_){
			wellName=wellName_;
			reads=reads_;
			tagName=tagName_;
		}
		
		@Override
		public String toString(){
			return wellName+"\t"+tagName+"\t"+reads;
		}
		
		/** Source well name where contamination originated */
		String wellName;
		/** Number of contaminating reads with this tag */
		long reads;
		/** Kapa tag identifier */
		String tagName;
		
	}
	
	//Ugly because it was retrofitted to support unsorted arrays and plate names
	/**
	 * Aggregates contamination statistics for a specific Kapa tag across multiple observations.
	 * Collects parts-per-million measurements from different plates and provides
	 * statistical analysis methods for contamination reporting.
	 */
	class TagData implements Comparable<TagData> {
		
		/**
		 * Constructs TagData for a specific tag and well combination.
		 * @param name_ Kapa tag identifier
		 * @param wellName_ Associated well name
		 */
		TagData(String name_, String wellName_){
			name=name_;
			wellName=wellName_;
		}
		
		/**
		 * Adds a contamination measurement for a specific tag and plate.
		 * Records both the parts-per-million value and associated plate name
		 * for later statistical analysis.
		 *
		 * @param tag Target tag that was contaminated
		 * @param ppmk Parts per million of contamination
		 * @param plate Plate name where measurement was taken
		 */
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
		
		/**
		 * Gets array of parts-per-million values for a specific tag.
		 * @param key Tag identifier to retrieve values for
		 * @param sort Whether to sort values in ascending order
		 * @return Array of PPM contamination values
		 */
		double[] getPpmArray(String key, boolean sort){
			ArrayList<Double> list=ppmMap.get(key);
			return toPpmArray(list, sort);
		}
		
		/**
		 * Gets array of plate names associated with contamination measurements.
		 * @param key Tag identifier to retrieve plate names for
		 * @param sort Whether to sort plate names alphabetically
		 * @return Array of plate names where contamination was observed
		 */
		String[] getPlateNameArray(String key, boolean sort){
			ArrayList<String> list=plateNameMap.get(key);
			return toPlateArray(list, sort);
		}
		
		/**
		 * Converts ArrayList of PPM values to fixed-size array.
		 * Pads with zeros to match expected observation count.
		 *
		 * @param list ArrayList of PPM values
		 * @param sort Whether to sort values in ascending order
		 * @return Fixed-size array of PPM values
		 */
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
		
		/**
		 * Converts ArrayList of plate names to fixed-size array.
		 * Pads with nulls to match expected observation count.
		 *
		 * @param list ArrayList of plate names
		 * @param sort Whether to sort names alphabetically
		 * @return Fixed-size array of plate names
		 */
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
		
		/** Kapa tag identifier */
		final String name;
		/** Associated well name for this tag */
		final String wellName;
		/** Number of times this tag has been observed across plates */
		int timesSeen=0;

		/** Maps tag names to lists of PPM contamination values */
		LinkedHashMap<String, ArrayList<Double>> ppmMap=new LinkedHashMap<String, ArrayList<Double>>(203);
		/** Maps tag names to lists of associated plate names */
		LinkedHashMap<String, ArrayList<String>> plateNameMap=new LinkedHashMap<String, ArrayList<String>>(203);
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file path containing plate names and lot numbers */
	private String in1=null;
	/** Output file path for contamination results */
	private String out1=null;
	
	/** URL prefix for JGI API calls */
	private String addressPrefix="https://rqc.jgi.doe.gov/api/plate_ui/page/";
	/** URL suffix for JGI API Kapa data requests */
	private String addressSuffix="/kapa spikein";//"/kapa spikein";
	
	/** Whether to output raw data instead of statistical summary */
	private boolean printRaw=false;
	
	/*--------------------------------------------------------------*/
	
	/** Number of input lines processed */
	private long linesProcessed=0;
	/** Number of output lines written */
	private long linesOut=0;
	/** Number of input bytes processed */
	private long bytesProcessed=0;
	/** Number of output bytes written */
	private long bytesOut=0;
	
	/** Maximum number of input lines to process */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file format specification */
	private final FileFormat ffin1;
	/** Output file format specification */
	private final FileFormat ffout1;
	
	/** JSON parser for processing API responses */
	private final JsonParser jp=new JsonParser();

	/** Maps Kapa tag names to aggregated contamination data */
	private final LinkedHashMap<String, TagData> tagMap=new LinkedHashMap<String, TagData>(203);
	/** Maps plate names to Plate objects for analysis */
	private final LinkedHashMap<String, Plate> plateMap=new LinkedHashMap<String, Plate>(203);
	
	/** Dummy well object used when source well is not found */
	final Well dummy=new Well("X", new JsonObject(), "X");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and errors */
	private PrintStream outstream=System.err;
	/** Whether to output detailed processing information */
	public static boolean verbose=false;
	/** Whether an error has occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
