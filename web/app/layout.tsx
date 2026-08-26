import React from 'react'
import './globals.css'

export const metadata = {
  title: 'Parallel Pangenome Graph Explorer',
  description: 'Comparing monolithic and parallel/reassembled pangenome graph construction. BCM SV Hackathon 2026.',
  openGraph: {
    title: 'Parallel Pangenome Graph Explorer',
    description: 'Can regional pangenome graphs be built in parallel and stitched back together?',
    type: 'website',
  },
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  )
}
