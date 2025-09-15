package shared;

public class Colors {
	
	public static void main(String[] args){
		System.out.println(format("Red", RED, false));
		System.out.println(format("Blue", BLUE, false));
		System.out.println(format("Green", GREEN, true));
		System.out.println(format("Yellow", YELLOW, true));
	}
	
	public static String format(String s, String color, boolean underline){
//		return color+(bold ? "[1m":"")+(bold ? "[4m":"")+s+RESET;
		return color+(underline ? UNDERLINE : "")+s+RESET;
	}
	
	public static byte[] format(byte[] s, String color, boolean underline){
		return format(new String(s), color, underline).getBytes();
	}
	
	public static String[] makeColorArray(){
		return new String[] {
				BRIGHT_RED, BRIGHT_GREEN, BRIGHT_YELLOW, BRIGHT_BLUE, BRIGHT_PURPLE, BRIGHT_CYAN, RED, GREEN, YELLOW, PURPLE, CYAN
		};
	}
	
	public static String[] makeDarkArray(){
		return new String[] {
				RED, GREEN, YELLOW, BLUE, PURPLE, CYAN,
		};
	}
	
	public static String[] makeBrightArray(){
		return new String[] {
				BRIGHT_RED, BRIGHT_GREEN, BRIGHT_YELLOW, BRIGHT_BLUE, BRIGHT_PURPLE, BRIGHT_CYAN
		};
	}

	public static boolean windows = Shared.WINDOWS;
	public static boolean skip = windows;
	public static String esc = windows ? "<ESC>" : "\u001B"; //Windows only supports colors with Win10+

	public static String RESET = skip ? "" : esc+"[0m";
	public static String UNDERLINE = skip ? "" : esc+"[4m";
	public static String BOLD = skip ? "" : esc+"[1m";
	
	public static String BLACK = skip ? "" : esc+"[30m";
	public static String RED = skip ? "" : esc+"[31m";
	public static String GREEN = skip ? "" : esc+"[32m";
	public static String YELLOW = skip ? "" : esc+"[33m";
	public static String BLUE = skip ? "" : esc+"[34m";
	public static String PURPLE = skip ? "" : esc+"[35m";
	public static String CYAN = skip ? "" : esc+"[36m";
	public static String WHITE = skip ? "" : esc+"[37m";
	
	public static String BRIGHT_BLACK = skip ? "" : esc+"[30;1m";
	public static String BRIGHT_RED = skip ? "" : esc+"[31;1m";
	public static String BRIGHT_GREEN = skip ? "" : esc+"[32;1m";
	public static String BRIGHT_YELLOW = skip ? "" : esc+"[33;1m";
	public static String BRIGHT_BLUE = skip ? "" : esc+"[34;1m";
	public static String BRIGHT_PURPLE = skip ? "" : esc+"[35;1m";
	public static String BRIGHT_CYAN = skip ? "" : esc+"[36;1m";
	public static String BRIGHT_WHITE = skip ? "" : esc+"[37;1m";

	public static String[] colorArray=makeColorArray();
	public static String[] darkArray=makeDarkArray();
	public static String[] BrightArray=makeBrightArray();
	
}
