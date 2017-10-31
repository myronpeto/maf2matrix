package org.panda.resource.tcga;

import org.panda.resource.Download;
import org.panda.utility.FileUtil;
import org.panda.utility.StringUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.*;

/**
 * Downloads TCGA data and analysis files from Broad Firehose.
 *
 * @author Ozgun Babur
 */
public class BroadDownloader
{
	private static final String BROAD_URL_PREFIX = "http://gdac.broadinstitute.org/runs/";
	private static final String BROAD_DATA_URL_PREFIX = BROAD_URL_PREFIX + "stddata__";
	private static final String BROAD_ANALYSIS_URL_PREFIX = BROAD_URL_PREFIX + "analyses__";

	//http://gdac.broadinstitute.org/runs/stddata__2015_04_02/data/GBM/201504022015040200.0.0.tar.gz

	public static void downloadAll(String date, String dir) throws IOException
	{
		List<String> codes = getStudyCodes(date);
//		codes.retainAll(Arrays.asList("FPPP"));

		for (String code : codes)
		{
			download(date, dir, code);
		}
	}

	public static void download(String date, String dir, String code) throws IOException
	{
		Map<ResourceType, String> pack = detectURLs(date, code);

		for (ResourceType type : ResourceType.values())
		{
			if (!pack.containsKey(type))
			{
				System.out.println("No " + type + " for " + code);
			}
		}

		if (pack.containsKey(ResourceType.CNA)) downloadCopyNumber(dir, code, pack.get(ResourceType.CNA));
		if (pack.containsKey(ResourceType.EXPRESSION)) downloadExpression(dir, code, pack.get(ResourceType.EXPRESSION));
		if (pack.containsKey(ResourceType.MUTATIONS)) downloadMutations(dir, code, pack.get(ResourceType.MUTATIONS));
		if (pack.containsKey(ResourceType.MUTSIG)) downloadMutsigScores(dir, code, pack.get(ResourceType.MUTSIG));
		if (pack.containsKey(ResourceType.RPPA)) downloadRPPA(dir, code, pack.get(ResourceType.RPPA));
	}

	private static Map<ResourceType, String> detectURLs(String date, String code) throws IOException
	{
		String d = date.replaceAll("_", "");
		Map<ResourceType, String> pack = new HashMap<ResourceType, String>();
		detectURLs(BROAD_ANALYSIS_URL_PREFIX + date + "/data/" + code + "/" + d, pack);
		detectURLs(BROAD_DATA_URL_PREFIX     + date + "/data/" + code + "/" + d, pack);
		return pack;
	}

	private static void detectURLs(String url, Map<ResourceType, String> pack) throws IOException
	{
		Scanner sc = new Scanner(new URL(url).openStream());
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			for (ResourceType type : ResourceType.values())
			{
				for (String match : type.content)
				{
					if (!pack.containsKey(type))
					{
						String pointer = StringUtil.fetch(line, match, "\"", "\"");
						if (pointer != null)
						{
							pointer = url + "/" + pointer;
							pack.put(type, pointer);
						}
					}
					else break;
				}
			}
		}
	}

	public static boolean downloadCopyNumber(String dir, String code, String url) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		String outFile = directory + "copynumber.txt";
		if (new File(outFile).exists()) return true;
		new File(directory).mkdirs();

		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "all_thresholded.by_genes", directory + "copynumber.txt"))
			{
				FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "amp_genes.conf_99.txt", directory + "scores-amplified.txt");
				FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "del_genes.conf_99.txt", directory + "scores-deleted.txt");

				return new File(tempFile).delete();
			}
		}
		return false;
	}
	public static boolean downloadMutsigScores(String dir, String code, String url) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		String outFile = directory + "scores-mutsig.txt";
		if (new File(outFile).exists()) return true;
		new File(directory).mkdirs();

		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "sig_genes.txt", outFile))
			{
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	public static boolean downloadExpression(String dir, String code, String url) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		String outFile = directory + "expression.txt";
		if (new File(outFile).exists()) return true;
		new File(directory).mkdirs();

		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "data.txt", outFile))
			{
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	public static boolean downloadMutations(String dir, String code, String url) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		String outFile = directory + "mutation.maf";
		if (new File(outFile).exists()) return true;
		new File(directory).mkdirs();

		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractAllEntriesContainingNameInTARGZFile(tempFile, "maf.txt", directory))
			{
				uniteMutations(directory, outFile);
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	private static void uniteMutations(String dir, String outFile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		File d = new File(dir);
		for (File file : d.listFiles())
		{
			if (file.getName().endsWith(".maf.txt"))
			{
				Scanner sc = new Scanner(file);
				while (sc.hasNextLine())
				{
					writer.write(sc.nextLine() + "\n");
				}
				sc.close();
				file.delete();
			}
		}
		writer.close();
	}

	public static boolean downloadRPPA(String dir, String code, String url) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		String outFile = directory + "rppa.txt";
		if (new File(outFile).exists()) return true;
		new File(directory).mkdirs();

		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, ".rppa.txt", outFile))
			{
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	public static List<String> getStudyCodes(String date) throws IOException
	{
		List<String> codes = new ArrayList<String>();
		Scanner sc = new Scanner(new URL(BROAD_DATA_URL_PREFIX + date + "/ingested_data.tsv").openStream());
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("Totals")) break;

			codes.add(line.substring(0, line.indexOf("\t")));
		}
		return codes;
	}

	enum ResourceType
	{
		MUTATIONS(".Mutation_Packager_Calls.Level_3."),
		CNA(".CopyNumber_Gistic2.Level_4."),
		EXPRESSION(".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3."),
//			       ".Merge_rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.Level_3."),
		MUTSIG(".MutSigNozzleReport2CV.Level_4."),
		RPPA(".RPPA_AnnotateWithGene.Level_3.");

		String[] content;

		ResourceType(String... content)
		{
			this.content = content;
		}
	}

	public static void main(String[] args) throws IOException
	{
//		List<String> codes = getStudyCodes("2015_08_21");
//		System.out.println(codes);

		download("2016_01_28", "/home/babur/Documents/TCGA/UVM", "UVM");

//		downloadAll("2016_01_28", "/home/babur/Documents/TCGA");
	}
}
