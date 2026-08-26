import { NextResponse } from "next/server"
import fs from "fs"
import path from "path"

export async function GET() {
  try {
    const p = path.join(process.cwd(), "public", "data", "latest.json")
    return NextResponse.json(JSON.parse(fs.readFileSync(p, "utf-8")))
  } catch {
    return NextResponse.json({
      data_mode: "synthetic",
      run: { run_id: "placeholder", pipeline_version: "0.1.0", mode: "placeholder" },
      metrics: {}, samples: [],
      message: "No data yet. Run: make demo",
    })
  }
}
