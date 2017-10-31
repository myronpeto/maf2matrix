package org.panda.resource.tcga;

import org.apache.commons.io.FileUtils;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.OrderedAveragePlot;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/9/2016.
 */
public class AlterationMatrixSeparator
{
	public static void main(String[] args) throws IOException
	{
//		printHistogram("/home/babur/Documents/mutex/BMEG/HNSC/TCGA-HNSC.tsv", 100);

		String base = "/home/babur/Documents/mutex/BMEG/HNSC/";
//		separate(base, "whole", new String[]{"0-50", "0-100", "0-200", "0-400"});

//		String dir = "/home/babur/Documents/mutex/TCGA/";
//		for (File file : new File(dir).listFiles())
//		{
//			System.out.println("file.getName() = " + file.getName());
//			if (file.getName().endsWith("keep")) continue;
////			separateNonOutliers(file.getPath(), "whole");
//			separateMutations(file.getPath(), "outliers-excluded");
//		}

//		for (File dir : new File("/home/babur/Documents/mutex/TCGA").listFiles())
//		{
//			if (!dir.isDirectory() || dir.getName().equals("ACC")) continue;
//
//			copySubdirs(new File("/home/babur/Documents/mutex/TCGA/ACC/whole"),
//				new File(dir.getPath() + "/whole"));
//		}

//		separate(base, base + "/whole", 2);
		separateNonOutliers(base, base + "/whole");

//		separateWithGeneAlterationStatus("/home/babur/Documents/mutex/TCGA/PanCan/1/1/DataMatrix.txt",
//			"/home/babur/Documents/mutex/TCGA/PanCan-TP53/intact/1/1/DataMatrix.txt", "TP53",
//			alt -> alt == 0);

//		copyNewAnalysisDir();
	}

	/**
	 *
	 * @param baseDir
	 * @param wholeDir
	 * @param chunks each dir name has to be in the format "a-b" where a and b are minimum and
	 *               maximum sample alteration counts, such like "20-70".
	 */
	public static void separate(String baseDir, String wholeDir, String[] chunks) throws IOException
	{
		for (String chunk : chunks)
		{
			int min = Integer.parseInt(chunk.substring(0, chunk.indexOf("-")));
			int max = Integer.parseInt(chunk.substring(chunk.indexOf("-") + 1));

			String chunkDir = baseDir + "/" + chunk;
			if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

			writeSubset(wholeDir + "/DataMatrix.txt", chunkDir + "/DataMatrix.txt",
				min, max);
		}
	}

	public static void separate(String baseDir, String wholeDir, int pieces) throws IOException
	{
		String inFile = wholeDir + "/DataMatrix.txt";
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);

		List<String> samples = new ArrayList<>(cnt.keySet());
		Collections.sort(samples, (o1, o2) -> cnt.get(o1).compareTo(cnt.get(o2)));

		for (int i = 1; i <= pieces; i++)
		{
			List<String> sub = samples.subList(
				(int) (samples.size() * (i - 1) / (double) pieces),
				((int) (samples.size() * i / (double) pieces)) - 1);

			String chunkDir = baseDir + "/" + i;
			if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

			writeSubset(inFile, chunkDir + "/DataMatrix.txt", new HashSet<>(sub));
		}
	}

	public static void separateNonOutliers(String baseDir, String wholeDir) throws IOException
	{
		String inFile = wholeDir + "/DataMatrix.txt";
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);

		List<Integer> list = new ArrayList<>(cnt.values());
		Collections.sort(list);
		int q1 = list.get((int) Math.floor(list.size() * 0.25));
		int q3 = list.get((int) Math.floor(list.size() * 0.75));

		double iqr = q3 - q1;

		double thr = q3 + (1.5 * iqr);

		Set<String> subset = new HashSet<>();

		for (String s : cnt.keySet())
		{
			if (cnt.get(s) <= thr) subset.add(s);
		}

		String chunkDir = baseDir + "/outliers-excluded";
		if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

		writeSubset(inFile, chunkDir + "/DataMatrix.txt", subset);
	}

	public static void separateMutations(String baseDir, String useDir) throws IOException
	{
		String inFile = baseDir + "/" + useDir + "/DataMatrix.txt";

		String newDir = baseDir + "/mutations-only";
		if (!new File(newDir).exists()) new File(newDir).mkdir();

		Scanner sc = new Scanner(new File(inFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(newDir + "/DataMatrix.txt"));
		writer.write(sc.nextLine());

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			line = line.replaceAll("\t2", "\t0").replaceAll("\t3", "\t0").
				        replaceAll("\t4", "\t1").replaceAll("\t5", "\t1");

			if (line.contains("\t1")) writer.write("\n" + line);
		}

		sc.close();
		writer.close();

		copySubdirs(new File(baseDir + "/" + useDir), new File(newDir));
	}

	public static void separateWithGeneAlterationStatus(String input, String output, String gene, AlterationSelector selector) throws IOException
	{
		writeSubset(input, output, selectSamplesWithGeneAlterationStatus(input, gene, selector));
	}

	/**
	 *
	 * @param inFile data matrix
	 * @param min minimum sample alteration count
	 * @param max maximum sample alteration count
	 */
	public static void writeSubset(String inFile, String outFile, int min, int max) throws IOException
	{
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);
		Set<String> samples = new HashSet<>();
		for (String s : cnt.keySet())
		{
			int c = cnt.get(s);

			if (c >= min && c <= max) samples.add(s);
		}

		writeSubset(inFile, outFile, samples);
	}

	public static void writeSubset(String inFile, String outFile, Set<String> samples) throws IOException
	{
		Scanner sc = new Scanner(new File(inFile));
		String line = sc.nextLine();
		String[] header = line.split("\t");

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		writer.write(header[0]);

		for (int i = 1; i < header.length; i++)
		{
			if (samples.contains(header[i])) writer.write("\t" + header[i]);
		}

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			writer.write("\n" + token[0]);
			for (int i = 1; i < token.length; i++)
			{
				if (samples.contains(header[i])) writer.write("\t" + token[i]);
			}
		}
		writer.close();

		copySubdirs(new File(inFile).getParentFile(), new File(outFile).getParentFile());
	}

	private static void copySubdirs(File dirFrom, File dirTo) throws IOException
	{
		for (File dir : dirFrom.listFiles())
		{
			if (dir.isDirectory())
			{
				FileUtils.copyDirectory(dir, new File(dirTo.getPath() + "/" + dir.getName()));
			}
		}
	}

	private static void copyNewAnalysisDir() throws IOException
	{
		String base = "/home/babur/Documents/mutex/TCGA/PanCan-TP53/intact/";
		File dirFrom = new File(base + "1/1");

		for (int i = 2; i < 10; i++)
		{
			for (int j = 1; j <= i; j++)
			{
				copySubdirs(dirFrom, new File(base + i + "/" + j));
			}
		}
	}

	public static void printHistogram(String file, int binSize) throws FileNotFoundException
	{
		Map<String, Integer> counts = readSampleAlterationCounts(file);
		Histogram h = new Histogram(binSize);
		h.setBorderAtZero(true);
		for (Integer i : counts.values())
		{
			h.count(i);
		}
		h.print();

		System.out.println();
		OrderedAveragePlot oap = new OrderedAveragePlot(counts.values(), true);
		oap.plot(30);
	}

	public static Map<String, Integer> readSampleAlterationCounts(String file) throws FileNotFoundException
	{
		Map<String, Integer> map = new HashMap<>();
		Scanner sc = new Scanner(new File(file));
		String line = sc.nextLine();
		line = line.substring(line.indexOf("\t") + 1);
		String[] samples = line.split("\t");
		for (String sample : samples)
		{
			map.put(sample, 0);
		}
		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			line = line.substring(line.indexOf("\t") + 1);
			String[] row = line.split("\t");
			for (int i = 0; i < row.length; i++)
			{
				if (!row[i].equals("0")) map.put(samples[i], map.get(samples[i]) + 1);
			}
		}
		sc.close();
		return map;
	}

	public static Set<String> selectSamplesWithGeneAlterationStatus(String file, String gene, AlterationSelector selector)
		throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(file));
		String[] header = sc.nextLine().split("\t");
		String[] row = null;
		while (sc.hasNextLine() && row == null)
		{
			String line = sc.nextLine();
			if (line.startsWith(gene + "\t"))
			{
				row = line.split("\t");
			}
		}
		if (row == null) throw new RuntimeException("Cannot find the gene: " + gene);

		Set<String> samples = new HashSet<>();
		for (int i = 1; i < header.length; i++)
		{
			if (selector.select(Integer.parseInt(row[i]))) samples.add(header[i]);
		}
		return samples;
	}

	interface AlterationSelector
	{
		boolean select(int alt);
	}
}
